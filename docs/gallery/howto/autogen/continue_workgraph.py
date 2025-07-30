"""
Continue a finished workgraph
=============================
"""


# %%
# There are two ways to continue a process in ``WorkGraph``:
#
# 1. Rerun the workgraph with modified inputs, relying on AiiDA's caching mechanism to avoid redundant calculations.
#    Visit the AiiDA documentation section on `caching <https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/provenance/caching.html#controlling-caching>`_ to learn more.
# 2. Load and restart a finished workgraph, allowing you to rerun it with modified inputs or even new tasks.
#    This feature is only available via the advanced usage paradigms.
#    For example, see how to do it using the :ref:`context manager <advanced:context-manager:continue-workflow>` paradigm.
