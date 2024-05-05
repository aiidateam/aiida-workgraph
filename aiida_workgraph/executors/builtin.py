from aiida.engine import WorkChain
from aiida import orm
from aiida.engine.processes.workchains.workchain import WorkChainSpec


class GatherWorkChain(WorkChain):
    @classmethod
    def define(cls, spec: WorkChainSpec) -> None:
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

    def gather(self) -> None:
        datas = self.inputs.datas.values()
        uuids = [data.uuid for data in datas]
        # uuids = gather(uuids)
        self.out("result", orm.List(uuids).store())
