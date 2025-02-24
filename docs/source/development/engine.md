# Engine

## Handle WorkGraph inputs
Use `metadata.workgraph_data` to save the `wgdata`. To do this we need to serialize the data, however, some of the input AiiDA data node is not stored.


## Updatable attributes

In `on_create` method, the `node` store itself, and its attributes become immutable. However, for the state of the tasks, we want to keep them updatable. To achieve this, we have to define the attributes that can be updated in the `_updatable_attributes` list.
