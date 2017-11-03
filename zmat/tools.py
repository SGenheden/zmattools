"""
Tools to generate and manipulate z-matrices
"""
from __future__ import division, print_function, absolute_import

def main_zmat() :
    """
    Entry point for commands associated with ELBA.
    The available tools are listed from generate.py and modify.py
    """

    import inspect
    import sys
    import zmat.commands as commands

    def __pred(c) :
        return inspect.isclass(c) and c.__module__ == "zmat.commands" and \
                issubclass(c,commands.ZmatCommand) and \
                c.__name__ != "ZmatCommand"

    cmdclasses = {name.lower():c
        for name,c in inspect.getmembers(sys.modules["zmat.commands"],__pred)}

    commands = {}
    helplines = []
    for cmdclass in sorted(cmdclasses) :
        cmdstr = cmdclasses[cmdclass].command_name()
        commands[cmdstr] = cmdclass
        helplines.append("%-25s - %s"%(cmdstr,
                            cmdclasses[cmdclass].descr()))

    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1].lower() == "-h") :
        help = """
-------------------------------------------------------
zmat - tools to generate and manipulate z-matrices
-------------------------------------------------------

Available commands:
%s
"""
        print(help%"\n".join(helplines))
    elif len(sys.argv) > 1 and sys.argv[1] in commands :
            clsname = commands[sys.argv[1]]
            cmdclasses[clsname]().execute()
    else :
        print ("Unrecognizable command: use -h for help.")
