## Summary

<!-- What changes, and why. Link the todo.md item(s) this closes, if any. -->

## Validation

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --workspace --all-targets --all-features -- -D warnings`
- [ ] `cargo test --workspace --all-features`
- [ ] `cargo deny check` (if dependencies changed)
- [ ] `no_std` builds still pass (if a no_std crate changed)

## Physics and references

<!-- For changes to formulas or constants: the source (standard, paper, or
     analytic result) and the tolerance the test enforces. -->

## Checklist

- [ ] New public items have rustdoc
- [ ] `todo.md` is updated, and ticked items refer to code that exists
