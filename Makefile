CLFAGS=-std=c99 -Wall -Werror -Wextra -lm
CC=gcc
PROJ=proj2
TEST=test.sh

.PHONY: tests clean debug

$(PROJ): $(PROJ).c
	$(CC) $(CLFAGS) $(PROJ).c -o $(PROJ)

debug: $(PROJ).c
	$(CC) $(CLFAGS) -DDEBUG $(PROJ).c -o $(PROJ)

tests: $(PROJ) tests/$(TEST)
	chmod +x tests/$(TEST)
	./tests/$(TEST) --show
	./tests/$(TEST) ./$(PROJ)

clean:
	rm -f $(PROJ)
