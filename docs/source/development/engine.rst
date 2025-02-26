====================
Engine
====================


In ``aiida-workgraph``, the ``WorkGraphEngine`` class is responsible for executing workflows. Unlike ``WorkChain``, which requires predefined inputs, outputs, and an execution outline for each workflow, ``WorkGraphEngine`` does not have built-in workflow definitions. Instead, it relies on a dynamic input namespace, ``workgraph_data``, allowing users to provide workflow information at runtime. This approach enables greater flexibility in constructing and executing workflows.

Defining a Workflow
-------------------
Users can define a workflow using:

1. The ``WorkGraph`` object (by adding tasks and links programmatically).
2. A YAML file.
3. A GUI (planned for future development).

In this document, we illustrate the process using the ``WorkGraph`` object.

Execution Flow
--------------

When executing a ``WorkGraph``, the following steps occur:

1. **Prepare Inputs**: Convert the ``WorkGraph`` to a dictionary (raw input data) and pass it to the ``WorkGraphEngine``.
2. **Create Process Node**: A ``WorkGraphNode`` is instantiated to represent the workflow execution.
3. **Extract and Store Workflow Data**:
   - Workflow structure (tasks, links, inputs, outputs, etc.).
   - Task states.
   - Task processes.
   - Task actions.
4. **Execute Workflow**: Tasks are processed, with their states, processes, and actions updated as needed.

Updatable Attributes
--------------------

Once a ``WorkGraphNode`` is created (via the ``on_create`` method), it is stored, making its attributes immutable. However, certain attributes must remain updatable during execution to reflect the workflow's progress.

- **State**: Represents the current status of a task.
- **Action**: Specifies interactions such as ``pause`` or ``play``.
- **Process**: Stores the primary key (``pk``) of the process node, as a task may launch multiple processes.

To allow these updates, the ``_updatable_attributes`` list defines the fields that can be modified during execution. This list is managed by the ``Sealable`` class, ensuring that only specified attributes can be changed.
