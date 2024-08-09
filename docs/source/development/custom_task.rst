==================
Custom Task
==================


Suggested naming rule for entry points:

- Format: {package_name}.{short_name_for_class}. If the plugin name starts with 'aiida', such as 'aiida-xxx', you may omit 'aiida' for brevity, e.g., 'xxx.add'. Aim to maintain consistency with the name used when registering with the AiiDA registry.
- The name is case-insensitive, thus `xxx.add` is the same as `xxx.Add`.
- Naming Style: Use snake case. For instance, use 'xxx.test_add' for the class 'TestAddTask'. Typically, the suffix 'task' is unnecessary, thus avoid using 'xxx.test_add_task'.

The rule is used for Task, Socket and Property.
