Data Serialization
===========================


In AiiDA
-------------------------------

In AiiDA, data needs to be serialized for two main reasons:
- to store it in the database, thus keep data provenance.
- to store intermediate data in the checkpoint, to be able to restart calculations.


In WorkGraph
-------------------------------
In workGraph, all input data are pass into the `wg` namespace. The `wg` namespace is a dynamic namespace. Since any data type can be passed into this namespace, so we don't add any vlidation on the data type. However, in AiiDA, when we generate the node graph, only the AiiDA data type will be shown in the graph, on another word, only the AiiDA data type can be linked inside the graph. That doesn't mean we lose the data provenance, since WorkGraph itself will not generate any data, and only the `calcfunction`, `CalcJob` will generate new data, and the input data of this process are required to be AiiDA data type. So the data provenance is still kept. Thus, one can always find the data provenance by checking the input data of the `calcfunction`, `CalcJob`.

One may argue that, why don't we serialize all data in the `wg` namespace? We choose not to do this for the following reasons:
- WorkGraph also support non-AiiDA component as node, that is to say, we can use any Python function as a node in the graph. These functions don't need AiiDA data as input, and we don't want to restrict the input data type of these functions.
- For AiiDA components (e.g. `CalcJob`, `WorkChain`), they may already define the serialization method for some input port explicity, thus, we need to respect this serialization method.

However, it still good to check if all data inside the `wg` namespace are json-serializable, to make sure that the checkpoint can be saved and loaded correctly.

PythonJob
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A special case is the `PythonJob`, which is a `CalcJob` that runs a Python function on a remote computer. The input data of the function to be fun is not required to be AiiDA data type, and the users are not required to provide the AiiDA data type as input. When the WorkGraph launches the `PythonJob`, it will serialize all input data of the function. However, if user use non-json-serializable data as input, the checkpoint will fail. Thus, we need to serialize all input data of the function when initializing the WorkGraph process, that is to say we need to serialize the inside the `wgdata`.
