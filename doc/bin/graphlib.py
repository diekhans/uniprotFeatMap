import os.path as osp
import graphviz as gv  # pip install graphviz

def check_output_file_type(parser, outfile):
    if _get_output_format(outfile) not in gv.FORMATS:
        fmts = ", ".join(gv.FORMATS)
        parser.error(f"output file must end in one of the format extensions: {fmts}: '{outfile}'")

def _get_output_format(outfile):
    return osp.splitext(outfile)[1][1:]

def _get_output_base(outfile):
    return osp.splitext(outfile)[0]

class GraphBuilder:
    def __init__(self, title, outfile, page=None):
        self.outfile = outfile
        fontname = "helvetica"
        graph_attr = {"fontname": fontname,
                      "label": title,
                      "labelloc": "top",
                      "fontsize": "25",
                      "rankdir": 'TB',
                      'ratio': 'auto'}
        graph_attr['ranksep'] = '0.5'
        if page is not None:
            graph_attr["page"] = page
        node_attr = {"fontname": fontname}
        edge_attr = {"fontname": fontname}
        self.graph = gv.Digraph(
            name=title,
            format=_get_output_format(self.outfile),
            graph_attr=graph_attr,
            node_attr=node_attr,
            edge_attr=edge_attr)
        self.data_files = set()
        self.stack = []

    def ensure_data_file(self, name):
        "create an data file node, if it does not already exist"
        if name not in self.data_files:
            self.graph.node(name, label=name,
                            shape="cylinder",
                            style="filled",
                            fillcolor="cyan",
                            fontcolor="black")
            self.data_files.add(name)

    def ensure_ext_data_file(self, name):
        "create an external data file node, if it does not already exist"
        if name not in self.data_files:
            self.graph.node(name, label=name,
                            shape="cylinder",
                            style="filled",
                            fillcolor="yellow",
                            fontcolor="black")
            self.data_files.add(name)

    def add_task(self, name, *, extins=(), inputs=(), outputs=()):
        "create a task node, inputs and outputs are other node names"
        label = name.replace(" ", "\\n")
        self.graph.node(name, label=label,
                        shape="box",
                        style="filled",
                        fillcolor="maroon",
                        fontcolor="white")
        for n in extins:
            self.ensure_ext_data_file(n)
            self.graph.edge(n, name)
        for n in inputs:
            self.ensure_data_file(n)
            self.graph.edge(n, name)
        for n in outputs:
            self.ensure_data_file(n)
            self.graph.edge(name, n)

    def add_program(self, name, *, inputs=(), outputs=()):
        self.graph.node(name,
                        label=name,
                        shape="trapezium",
                        style="filled",
                        fillcolor="darkgreen",
                        fontcolor="white")
        for n in inputs:
            self.ensure_data_file(n)
            self.graph.edge(n, name)
        for n in outputs:
            self.ensure_data_file(n)
            self.graph.edge(name, n)

    def _push(self, new_graph):
        self.stack.append(self.graph)
        self.graph = new_graph

    def push_cluster(self, name):
        invis = {"style": "invis"}
        self._push(self.graph.subgraph(name="cluster_" + name, graph_attr=invis))

    def push_group(self, name):
        self._push(self.graph.subgraph(name=name))

    def pop(self):
        self.graph = self.stack.pop()

    def render(self, keep=False):
        self.graph.render(_get_output_base(self.outfile), cleanup=not keep)
