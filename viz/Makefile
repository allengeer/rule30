CC ?= cc
CFLAGS ?= -O3 -march=native -Wall -Wextra
BUILD = build

all: $(BUILD)/r30gen $(BUILD)/centre_prep $(BUILD)/libr30.dylib

$(BUILD)/r30gen: r30gen.c r30core.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ r30gen.c

$(BUILD)/centre_prep: centre_prep.c | $(BUILD)
	$(CC) $(CFLAGS) -o $@ centre_prep.c

$(BUILD)/libr30.dylib: r30lib.c r30core.h | $(BUILD)
	$(CC) $(CFLAGS) -shared -fPIC -o $@ r30lib.c

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)
