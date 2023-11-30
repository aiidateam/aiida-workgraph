from aiida.cmdline.utils.common import get_workchain_report
import aiida

aiida.load_profile()


def test_run_order(arithmetic_add):
    """Test the order.
    Nodes should run in parallel and only depend on the input nodes."""
    from aiida.orm import load_code, Int
    from aiida_worktree import WorkTree

    code = load_code("add@localhost")
    x = Int(2)
    wt = WorkTree(name="test_run_order")
    adds = []
    for i in range(6):
        temp = wt.nodes.new(arithmetic_add, f"add{i}", x=x, y=Int(i), code=code)
        if i == 0:
            temp.set({"metadata.options.sleep": 15})
        else:
            temp.set({"metadata.options.sleep": 1})
        adds.append(temp)
    wt.links.new(adds[0].outputs["sum"], adds[2].inputs["x"])
    wt.links.new(adds[1].outputs["sum"], adds[3].inputs["x"])
    wt.links.new(adds[3].outputs["sum"], adds[4].inputs["x"])
    wt.links.new(adds[2].outputs["sum"], adds[5].inputs["x"])
    wt.links.new(adds[4].outputs["sum"], adds[5].inputs["y"])
    wt.submit(wait=True)
    wt.process.pause()
    # get the report
    report = get_workchain_report(wt.process, "REPORT")
    lines = report.splitlines()
    steps = {}
    for i, line in enumerate(lines):
        if "Run node:" in line:
            print(line)
            name = line.split("Run node:")[1].split(",")[0].strip()
            steps[name] = i
    # first run add0, add1
    # then add3
    # then add4
    # then add2
    # finally add5
    assert steps["add2"] > steps["add4"]
