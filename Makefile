.PHONY: build clean edit install uninstall reinstall

build:
	dune build --profile=release @install -j 16

clean:
	rm -rf _build

edit:
	emacs src/*.ml TODO commands.sh &

install: build
	cp bin/fasmifra_fragment.py `opam config var bin`/
	dune install

uninstall:
	dune uninstall
	rm -f `opam config var bin`/fasmifra_fragment.py

reinstall: uninstall install
