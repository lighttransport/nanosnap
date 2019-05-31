rmdir /q /s build

cmake -G "Visual Studio 15 2017" -A x64 -Bbuild -H.
