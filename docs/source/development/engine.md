# Engine

In WorkGraph engine, the tasks and their inputs and outputs are passed as inputs of the WorkGraphEngine. This is different from the WorkChain, where the developers hardcode the input and output port, and the outline for every workflow.

When running a WorkGraph, there are the following basic steps:

- export the WorkGraph as a dictionary (raw inputs), and pass it as input to the WorkGraphEngine
- a `WorkGraphNode` is created by the process
- get the WorkGraph data from the raw inputs, and save it as attributes of the `WorkGraphNode`, this includes:
    - WorkGraph data (takss, links, inputs, outputs etc)
    - task states
    - task processes
    - task actions
- run the WorkGraph, and update the task states, processes, and actions


## Updatable attributes
When the `WorkGraphNode` is created (in the `on_create` method), it is also stored, and its attributes become immutable. However, some of the data need updating during the execution of the WorkGraph, this includes:

- state of the task
- action (e.g., pause, reset) of the task
- process of the task, because


To achieve this, we have to define the attributes that can be updated in the `_updatable_attributes` list.
