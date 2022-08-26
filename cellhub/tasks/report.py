'''
Report.py
=========

Functions and variable definitions to help with generation of latex reports

Code
====

'''

import types

# --------------------------------- Templates -------------------------------- #

template = {"figure": '''\
            \\begin{figure}[H]
            \\includegraphics[width=%(width)s\\textwidth,height=%(height)s\\textheight,keepaspectratio]{{{%(path)s}}}
            \\caption{%(caption)s}
            \\end{figure}
            ''',
            "section": '''\\section{%(title)s}''',
            "subsection": '''\\subsection{%(title)s}'''
}

template = types.SimpleNamespace(**template)