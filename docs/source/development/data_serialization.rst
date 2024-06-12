==================
Data Serialization
==================

In AiiDA
--------

Data serialization in AiiDA is critical for several reasons:

- **Storing data in the database**: This is important for maintaining data provenance.

- **Storing intermediate data in checkpoints**: This enables the restart of calculations.

In WorkGraph
------------

In WorkGraph, all input data are passed into the ``wg`` namespace. The ``wg`` namespace is dynamic and can accept any data type, without enforcing validation on the data types. However, in AiiDA, when generating the node graph, only AiiDA data types are displayed, meaning only AiiDA data types can be linked within the graph. This does not imply a loss of data provenance, as WorkGraph itself does not generate any data. Only ``calcfunction`` and ``CalcJob`` generate new data, and the input data for these processes must be AiiDA data types to preserve data provenance. Data provenance can always be traced by checking the input data of the ``calcfunction`` and ``CalcJob``.

There are reasons why we don't serialize all data in the ``wg`` namespace:

- **Flexibility for non-AiiDA components**: WorkGraph supports non-AiiDA components as nodes, meaning any Python function can be used as a node in the graph. These functions do not require AiiDA data as input, allowing for a variety of data types.

- **Respecting existing serialization methods**: For AiiDA components (e.g., ``CalcJob``, ``WorkChain``), some input ports may have explicitly defined serialization methods, which must be respected.

However, ensuring that all data within the ``wg`` namespace are JSON-serializable is beneficial to guarantee that checkpoints can be saved and loaded correctly.

PythonJob
---------

``PythonJob`` is a special case of ``CalcJob`` that runs a Python function on a remote computer. The input data for the function does not need to be of AiiDA data type, and users are not required to provide AiiDA data types as input. When WorkGraph launches the ``PythonJob``, it serializes all input data for the function. However, if users provide non-JSON-serializable data as input, the checkpoint will fail. Thus, it is necessary to serialize all input data of the function when initializing the WorkGraph process.
