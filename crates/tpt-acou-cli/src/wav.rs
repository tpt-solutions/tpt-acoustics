//! Minimal WAV reader for the CLI: 16-bit PCM and 32-bit IEEE float, any
//! channel count (the first channel is returned).

/// Samples and sample rate of a WAV file, first channel only.
#[derive(Debug, Clone, PartialEq)]
pub struct WavData {
    /// Samples, nominally in [-1, 1].
    pub samples: Vec<f64>,
    /// Sample rate in hertz.
    pub sample_rate: u32,
}

fn u16_at(bytes: &[u8], at: usize) -> Option<u16> {
    Some(u16::from_le_bytes([*bytes.get(at)?, *bytes.get(at + 1)?]))
}

fn u32_at(bytes: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes([
        *bytes.get(at)?,
        *bytes.get(at + 1)?,
        *bytes.get(at + 2)?,
        *bytes.get(at + 3)?,
    ]))
}

/// Parses a WAV file.
///
/// # Errors
/// Returns a message for a malformed RIFF structure, a missing `fmt ` or
/// `data` chunk, or an unsupported encoding (only PCM16 and float32 are read).
pub fn parse_wav(bytes: &[u8]) -> Result<WavData, String> {
    if bytes.len() < 12 || &bytes[0..4] != b"RIFF" || &bytes[8..12] != b"WAVE" {
        return Err("not a RIFF/WAVE file".into());
    }
    let mut pos = 12;
    let mut format: Option<(u16, u16, u32, u16)> = None; // (tag, channels, rate, bits)
    let mut data: Option<&[u8]> = None;
    while pos + 8 <= bytes.len() {
        let id = &bytes[pos..pos + 4];
        let size = u32_at(bytes, pos + 4).ok_or("truncated chunk header")? as usize;
        let body_start = pos + 8;
        let body_end = body_start.checked_add(size).ok_or("chunk size overflows")?;
        if body_end > bytes.len() {
            return Err("chunk extends past end of file".into());
        }
        let body = &bytes[body_start..body_end];
        if id == b"fmt " {
            let tag = u16_at(body, 0).ok_or("short fmt chunk")?;
            let channels = u16_at(body, 2).ok_or("short fmt chunk")?;
            let rate = u32_at(body, 4).ok_or("short fmt chunk")?;
            let bits = u16_at(body, 14).ok_or("short fmt chunk")?;
            format = Some((tag, channels, rate, bits));
        } else if id == b"data" {
            data = Some(body);
        }
        // Chunks are word-aligned.
        pos = body_end + (size & 1);
    }
    let (tag, channels, rate, bits) = format.ok_or("missing fmt chunk")?;
    let data = data.ok_or("missing data chunk")?;
    if channels == 0 {
        return Err("zero channels".into());
    }
    let stride = usize::from(channels);
    let samples: Vec<f64> = match (tag, bits) {
        (1, 16) => data
            .chunks_exact(2 * stride)
            .map(|frame| f64::from(i16::from_le_bytes([frame[0], frame[1]])) / 32768.0)
            .collect(),
        (3, 32) => data
            .chunks_exact(4 * stride)
            .map(|frame| f64::from(f32::from_le_bytes([frame[0], frame[1], frame[2], frame[3]])))
            .collect(),
        _ => return Err(format!("unsupported WAV encoding: tag {tag}, {bits} bits")),
    };
    Ok(WavData {
        samples,
        sample_rate: rate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_acoustics::ray::rir::write_wav_f32;

    #[test]
    fn reads_back_a_float_wav_written_by_the_library() {
        let original = [0.0, 0.5, -0.25, 1.0];
        let bytes = write_wav_f32(&original, 48_000);
        let back = parse_wav(&bytes).unwrap();
        assert_eq!(back.sample_rate, 48_000);
        for (a, b) in back.samples.iter().zip(original) {
            assert!((a - b).abs() < 1e-7, "{a} vs {b}");
        }
    }

    #[test]
    fn reads_16_bit_pcm_and_takes_the_first_channel() {
        // Stereo, 2 frames: left = 16384 and -16384, right = 0.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"RIFF");
        bytes.extend_from_slice(&(36 + 8u32).to_le_bytes());
        bytes.extend_from_slice(b"WAVE");
        bytes.extend_from_slice(b"fmt ");
        bytes.extend_from_slice(&16u32.to_le_bytes());
        bytes.extend_from_slice(&1u16.to_le_bytes());
        bytes.extend_from_slice(&2u16.to_le_bytes());
        bytes.extend_from_slice(&44_100u32.to_le_bytes());
        bytes.extend_from_slice(&(44_100u32 * 4).to_le_bytes());
        bytes.extend_from_slice(&4u16.to_le_bytes());
        bytes.extend_from_slice(&16u16.to_le_bytes());
        bytes.extend_from_slice(b"data");
        bytes.extend_from_slice(&8u32.to_le_bytes());
        for v in [16384i16, 0, -16384, 0] {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        let wav = parse_wav(&bytes).unwrap();
        assert_eq!(wav.sample_rate, 44_100);
        assert_eq!(wav.samples.len(), 2);
        assert!((wav.samples[0] - 0.5).abs() < 1e-9);
        assert!((wav.samples[1] + 0.5).abs() < 1e-9);
    }

    #[test]
    fn rejects_non_wav_and_truncated_files() {
        assert!(parse_wav(b"hello").is_err());
        let mut bytes = write_wav_f32(&[0.1, 0.2], 8_000);
        bytes.truncate(bytes.len() - 3);
        assert!(parse_wav(&bytes).is_err());
    }
}
