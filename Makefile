CLFAGS=-std=c99 -Wall -Werror -Wextra -lm -DDEBUG
CC=gcc
PROJ=proj2
TEST=tests.sh
TEST_OUTPUT=output.txt
TEST_PATTERN=pattern.txt

.PHONY: tests clean

$(PROJ): $(PROJ).c
	$(CC) $(CLFAGS) $(PROJ).c -o $(PROJ)

tests/$(TEST_OUTPUT): tests/$(TEST) tests/$(TEST_PATTERN)
	rm -f tests/$(TEST_OUTPUT)
	chmod +x tests/$(TEST)
	cd tests && ./$(TEST)

tests: $(PROJ) tests/$(TEST_OUTPUT)
	diff -u tests/$(TEST_OUTPUT) tests/$(TEST_PATTERN)

clean:
	rm -f $(PROJ)
	rm -f tests/output.txt
