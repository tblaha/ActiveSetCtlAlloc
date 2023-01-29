
tests:
	cd ./harness && make all

verify: tests
	./harness/main.o verify 0
	./harness/main.o verify 1
	./harness/main.o verify 2
	./harness/main.o verify 3

clean:
	cd ./harness && make clean
