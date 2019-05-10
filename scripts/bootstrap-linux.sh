rm -rf build

cmake -Bbuild -H. -DSANITIZE_ADDRESS=On
