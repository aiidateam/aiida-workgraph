from aiida.engine import WorkChain, calcfunction
from aiida import orm


@calcfunction
def gather(inputs):
    """Add node."""
    return inputs.clone()


class GatherWorkChain(WorkChain):
    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        spec.input_namespace(
            "datas",
            dynamic=True,
            help=('Dynamic namespace for the datas, "{key}" : {Data}".'),
        )
        spec.outline(
            cls.gather,
        )
        spec.output(
            "result",
            valid_type=orm.List,
            required=True,
            help="A list of the uuid of the outputs.",
        )

    def gather(self):
        datas = self.inputs.datas.values()
        uuids = [data.uuid for data in datas]
        # uuids = gather(uuids)
        self.out("result", orm.List(uuids).store())


if __name__ == "__main__":
    from aiida.orm import Int
    from aiida.engine import run
    from aiida import load_profile

    load_profile()
    run(GatherWorkChain, datas={"a": Int(1), "b": Int(2)})
