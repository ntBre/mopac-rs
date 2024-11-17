clippy:
	cargo clippy --workspace --tests

test:
	cargo test --workspace --no-fail-fast -- --test-threads=1

doc:
	cargo doc --open
