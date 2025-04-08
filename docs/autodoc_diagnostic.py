"""
Sphinx extension that acts as a autodoc patch for documenting Indicator instances.

By default, diagnostic instances are skipped by autodoc because their subclass is not a builtin type of python.

Based on https://github.com/powerline/powerline/blob/83d855d3d73498c47553afeba212415990d95c54/docs/source/powerline_autodoc.py
"""

from __future__ import annotations

from sphinx.domains.python import PyFunction, PyXRefRole
from sphinx.ext import autodoc

from valenspy.diagnostic import Diagnostic


class DiagnosticDocumenter(autodoc.FunctionDocumenter):
    objtype = "diagnostic"

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        return isinstance(member, Diagnostic) or issubclass(type(member), Diagnostic)

    def get_doc(self):
        return super().get_doc()


class DiagnosticDirective(PyFunction):
    pass


def setup(app):
    app.add_autodocumenter(DiagnosticDocumenter)
    app.add_directive_to_domain("py", "diagnostic", DiagnosticDirective)
    app.add_role_to_domain("py", "diagnostic", PyXRefRole())