# numerical_mathematics_methods

[![License](https://img.shields.io/hexpm/l/plug.svg)](https://github.com/Vasar007/numerical_mathematics_methods/blob/master/LICENSE)
[![Build Status](https://travis-ci.com/Vasar007/numerical_mathematics_methods.svg?branch=master)](https://travis-ci.com/Vasar007/numerical_mathematics_methods)
[![CodeFactor](https://www.codefactor.io/repository/github/vasar007/numerical_mathematics_methods/badge)](https://www.codefactor.io/repository/github/vasar007/numerical_mathematics_methods)

Implementation of some numerical mathematics methods.

## Compiling

This project is compiled by Clang v5.0.0 and parameters: *clang++ -std=c++1z -O1 -Wall -Wextra -pedantic -Xclang -flto-visibility-public-std -fconstexpr-steps=1000000000*.

Also it would be compiled by GCC v7.1 and higher.

MSVC from the v19.14 (VS v15.7) has some issues with constexpr statements and sometimes "Internal compiler error" occurs. But code could be compiled too.

## License information

This project is licensed under the terms of the [Apache License 2.0](LICENSE).
