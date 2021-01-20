# About

The GenUI framework provides both API and GUI for molecular generators, QSAR modelling and chemical space visualization. This repository contains the implementation of the backend application used in the GenUI platform. 

Molecular generators can be integrated
as self-contained Python packages and their inputs, outputs and data storage can 
be managed through various REST API endpoints implemented in this 
repository. It is possible to upload bioactivity datasets and also download 
them from public sources. QSAR and generative models can be built from this 
data as well as 2D representations of chemical space. The GenUI [web frontend](https://github.com/martin-sicho/genui-gui) uses the REST
 services implemented here to allow simpler access to all these features. 
 You can find out more about how to extend 
and use the GenUI Python API in the [documentation](https://martin-sicho.github.io/genui/docs/index.html).

You can also directly deploy and test existing open source 
GenUI ecosystem as Docker images 
(more details [here](https://github.com/martin-sicho/genui-docker)).

# Contact

If you have any questions, problems or feature requests, do not hesitate to 
create a GitHub issue. This application is still very much in development so we 
appreciate any feedback. We are also interested in collaborations to integrate 
more molecular generators in the platform so do not hesitate to contact us
if you would like to contribute.

# License

This software is distributed under [MIT license](./src/LICENSE.md)