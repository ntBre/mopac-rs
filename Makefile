clippy:
	cargo clippy --workspace --tests

test:
	cargo test --workspace --no-fail-fast
