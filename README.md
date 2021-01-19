# About

This repository contains the implementation of the backend application used in the GenUI framework. The GenUI framework provides both API and GUI for molecular generators. Molecular generators can be integrated
as self-contained Python packages and their inputs, outputs and data storage can 
be managed through various REST API endpoints implemented in this 
repository. It is possible to upload bioactivity datasets and also download 
them from public sources. QSAR and generative models can be built from this 
data as well as 2D representations of chemical space. The GenUI [web frontend](https://github.com/martin-sicho/genui-gui) uses the REST
 services implemented here to allow simple access to all these features
 to less experienced users as well and data visualization. 
 You can find out more about how to extend 
and use the GenUI Python API in the [documentation](https://martin-sicho.github.io/genui/docs/index.html).

You can also directly deploy and test existing open source 
GenUI ecosystem as Docker images 
(more details [here](https://github.com/martin-sicho/genui-docker)).

# Contact

If you have any questions, problems or feature requests, do not hesitate to 
create a GitHub issue. This application is still very much in development so we 
appreciate any feedback. 

# License

This software is distributed under the MIT license with the following terms:

```
Copyright 2021 Martin Šícho

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```