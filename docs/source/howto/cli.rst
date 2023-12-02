.. _topics:cli:

**********************
Command line interface
**********************
As a supplement command for `verdi`, the command line interface utility for AiiDA-WorkTree is called ``worktree``.
This section explains the basic concepts that apply to all ``worktree`` commands.

.. _topics:cli:parameters:

Parameters
==========
Parameters to ``worktree`` commands come in two flavors:

* Arguments: positional parameters, e.g. ``123`` in ``worktree tree show 123``
* Options: announced by a flag (e.g. ``-f`` or ``--flag``), potentially followed by a value. E.g. ``worktree tree list --limit 10`` or ``worktree tree -h``.


.. _topics:cli:help_strings:

Help strings
============
Append the ``--help`` option to any verdi (sub-)command to get help on how to use it.
For example, ``worktree tree kill --help`` shows::

    Usage: worktree tree kill [OPTIONS] [PROCESSES]...

        Kill running processes.

    Options:
        -t, --timeout FLOAT  Time in seconds to wait for a response before timing
                             out.  [default: 5.0]
        --wait / --no-wait   Wait for the action to be completed otherwise return as
                             soon as it's scheduled.
        -h, --help           Show this message and exit.
