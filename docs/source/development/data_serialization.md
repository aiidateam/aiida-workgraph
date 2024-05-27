Data Serialization
===================

In AiiDA
--------

In AiiDA, data serialization is crucial for two main reasons:
- To store data in the database, thereby maintaining data provenance.
- To store intermediate data in checkpoints, enabling the restart of calculations.

In WorkGraph
------------

In WorkGraph, all input data are passed into the `wg` namespace. The `wg` namespace is dynamic and can accept any data type, so we do not enforce validation on the data types. However, in AiiDA, when generating the node graph, only AiiDA data types are displayed. In other words, only AiiDA data types can be linked within the graph. This does not mean we lose data provenance, as WorkGraph itself does not generate any data. Only `calcfunction` and `CalcJob` generate new data, and the input data for these processes are required to be AiiDA data types. Therefore, data provenance is preserved. One can always trace data provenance by checking the input data of the `calcfunction` and `CalcJob`.

One may question why we don't serialize all data in the `wg` namespace. We choose not to for the following reasons:
- WorkGraph supports non-AiiDA components as nodes, meaning any Python function can be used as a node in the graph. These functions do not require AiiDA data as input, and we do not want to restrict their input data types.
- For AiiDA components (e.g., `CalcJob`, `WorkChain`), serialization methods for some input ports may already be explicitly defined, and we need to respect these serialization methods.

However, it is still beneficial to ensure that all data inside the `wg` namespace are JSON-serializable to guarantee that checkpoints can be saved and loaded correctly.

PythonJob
---------

A special case is the `PythonJob`, which is a `CalcJob` that runs a Python function on a remote computer. The input data for the function does not need to be of AiiDA data type, and users are not required to provide AiiDA data types as input. When WorkGraph launches the `PythonJob`, it serializes all input data for the function. However, if users provide non-JSON-serializable data as input, the checkpoint will fail. Thus, we need to serialize all input data of the function when initializing the WorkGraph process.
