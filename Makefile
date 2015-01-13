
VERSION=1.0

TARGET_BIN=/bin
SRC=./bin

install:
        cp -v ${SRC}/* ${TARGET_BIN}

clean:
	  rm /bin/ngseas* && rm /bin/


