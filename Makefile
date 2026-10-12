# Binary name
BINARY_DIR=golooper

# Directories
TEST_OUTPUT_DIR=test_output
OUT_DIR=out

# Build target
.PHONY: build
build: clean
	go build -o golooper

# Clean target
.PHONY: clean
clean:
		rm -f $(BINARY_DIR)
	rm -rf $(TEST_OUTPUT_DIR)
	rm -rf $(OUT_DIR)

# Test target
.PHONY: test
test: build
	# Run unit tests
	go test ./...

	# Create test output directory
	mkdir -p $(TEST_OUTPUT_DIR)

	# Run test execution and save output
	./$(BINARY_DIR) ensemble -f res/pfc53_full.fa -o $(OUT_DIR)/results --localaverageenergy > $(TEST_OUTPUT_DIR)/test_execution.log 2>&1
