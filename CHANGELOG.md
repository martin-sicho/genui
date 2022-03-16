# Change Log

Previous version: `0.0.0.alpha1`

Current version: `0.0.0.alpha2`

## Changes

- updated all dependencies (Django 4, Celery 5, RDKit 2021.03, ...) to improve security and get the latest features (closes #8, closes #7)
- added an option to create exports for compound sets, supported file types:
    - SDF Files
- integrated the [engstmaps](./src/genui/maps/extensions/engstmaps) package with more chemical space projections

## Fixes
- pending tasks without a result are now showed correctly with others
- fixed potential database key errors when the `getDjangoModel` methods of the discovery API are called
