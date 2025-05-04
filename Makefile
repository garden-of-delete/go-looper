# Binary name
BINARY_NAME=golooper

# Directories
TEST_OUTPUT_DIR=test_output
OUT_DIR=out

# Build target
.PHONY: build
build: clean
	go build -o $(BINARY_NAME)

# Clean target
.PHONY: clean
clean:
	rm -f $(BINARY_NAME)
	rm -rf $(TEST_OUTPUT_DIR)

# Test target
.PHONY: test
test: build
	# Run unit tests
	go test ./...
	
	# Create test output directory
	mkdir -p $(TEST_OUTPUT_DIR)
	
	# Run test execution and save output
	./$(BINARY_NAME) perloop -f res/pfc53_full.fa -o $(OUT_DIR)/results -m 3 -s 0.07 -C > $(TEST_OUTPUT_DIR)/test_execution.log 2>&1

# Default target
.PHONY: all
all: test git