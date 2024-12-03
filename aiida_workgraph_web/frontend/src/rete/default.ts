import { NodeEditor, GetSchemes, ClassicPreset } from "rete";
import { AreaPlugin, AreaExtensions } from "rete-area-plugin";
import {
  ConnectionPlugin,
  Presets as ConnectionPresets
} from "rete-connection-plugin";
import { ReactPlugin, Presets, ReactArea2D } from "rete-react-plugin";
import { ScopesPlugin, Presets as ScopesPresets } from "rete-scopes-plugin";
import { MinimapExtra, MinimapPlugin } from "rete-minimap-plugin";
import {
  ContextMenuPlugin,
  Presets as ContextMenuPresets,
  ContextMenuExtra
} from "rete-context-menu-plugin";

import {
  AutoArrangePlugin,
  Presets as ArrangePresets,
  ArrangeAppliers
} from "rete-auto-arrange-plugin";
import AtomsItem from './AtomsItem.js'; // Adjust the path as necessary

// May move these interfaces to a separate file
interface NodeInput {
  name: string;
  identifier: string;
  properties: any;
  // Add other properties of NodeInput here if needed
}

interface NodeOutput {
  name: string;
  // Add other properties of NodeOutput here if needed
}

interface NodeData {
  label: string;
  inputs: NodeInput[];
  outputs: NodeOutput[];
}
interface LinkData {
  from_node: string;
  from_socket: string;
  to_node: string;
  to_socket: string;
  // Add other properties of LinkData here if needed
}

interface NodeMap {
  [key: string]: Node; // Assuming `Node` is the correct type for the nodes in nodeMap
}

class AtomsControl extends ClassicPreset.Control {
  constructor(public data: any) {
    super();
  }
}

export async function loadJSON(editor: NodeEditor<any>, area: any, workgraphData: any) {

  // Adding nodes based on workgraphData
  const nodeMap: NodeMap = {}; // To keep track of created nodes for linking
  for (const nodeId in workgraphData.nodes) {
    const nodeData = workgraphData.nodes[nodeId];
    const node = createDynamicNode(nodeData);
    await editor.addNode(node);
    nodeMap[nodeId] = node; // Storing reference to the node
  }
  // Adding connections based on workgraphData
  workgraphData.links.forEach(async (link: LinkData) => { // Specify the type of link here
    const fromNode = nodeMap[link.from_node];
    const toNode = nodeMap[link.to_node];
    if (fromNode && toNode) {
        await editor.addConnection(new Connection(fromNode, link.from_socket, toNode, link.to_socket));
    }
  });

  // Add while zones
  console.log("Adding while zone: ");
  for (const nodeId in workgraphData.nodes) {
    const nodeData = workgraphData.nodes[nodeId];
    const node_type = nodeData['node_type'];
    if (node_type === "WHILE" || node_type === "IF" || node_type === "ZONE") {
      // find the node
      const node = nodeMap[nodeData.label];
      const children = nodeData['children'];
      // find the id of all nodes in the editor that has a label in while_zone
      for (const nodeId in children) {
        const node1 = nodeMap[children[nodeId]];
        console.log("Setting parent of node", node1, "to", node);
        node1.parent = node.id;
        area.update('node', node1.id);
      }
      area.update('node', node.id);
    }
  }
}


class Node extends ClassicPreset.Node <
{ [key in string]: ClassicPreset.Socket },
{ [key in string]: ClassicPreset.Socket },
{
  [key in string]:
    | AtomsControl
    | ClassicPreset.Control
    | ClassicPreset.InputControl<"number">
    | ClassicPreset.InputControl<"text">;
}
> {
  width = 180;
  height = 100;
  parent?: string;
}
class Connection<N extends Node> extends ClassicPreset.Connection<N, N> {}

type Schemes = GetSchemes<Node, Connection<Node>>;
type AreaExtra = ReactArea2D<any> | MinimapExtra | ContextMenuExtra;


export function addControls(editor: NodeEditor<any>, area: any,
          workgraphData: any) {
  // resize the node based on the max length of the input/output names
  let maxSocketNameLength = 0;
  const nodes = editor.getNodes();
  // loop all node and add the input controls
  nodes.forEach((node) => {
    const nodeData = workgraphData.nodes[node.label];
    nodeData.inputs.forEach((input: NodeInput) => {
        const inp = node.inputs[input.name];
        // console.log("Adding control for input", input.name, "with identifier", input.identifier);
        if (input.identifier === "workgraph.int" || input.identifier === "workgraph.float" || input.identifier === "workgraph.aiida_int" || input.identifier === "workgraph.aiida_float") {
          inp.addControl(
            new ClassicPreset.InputControl("number", { initial: nodeData.properties[input.name].value, readonly: true })
          );
        }
        if (input.identifier === "workgraph.string" || input.identifier === "workgraph.aiida_string") {
          inp.addControl(
            new ClassicPreset.InputControl("text", { initial: nodeData.properties[input.name].value, readonly: true })
          );
        }
        if (input.identifier === "workgraph.aiida_structuredata") {
          inp.addControl(new AtomsControl(nodeData.properties[input.name].value));
          node.width += 0;
          node.height += 150;
    }
    area.update('node', node.id);
  })
  });
}

export function removeControls(editor: NodeEditor<any>, area: any, workgraphData: any) {
  const nodes = editor.getNodes();
  nodes.forEach((node) => {
      const nodeData = workgraphData.nodes[node.label];
      nodeData.inputs.forEach((input: NodeInput) => {
        const inp = node.inputs[input.name];
        inp.control = null;
        const identifier = input.identifier;
        if (identifier === "workgraph.aiida_structuredata") {
          node.width -= 0;
          node.height -= 150;
        }
    });
    area.update('node', node.id);
  });
}


function createDynamicNode(nodeData: any) {
  const node = new Node(nodeData.label);
  // resize the node based on the max length of the input/output names
  let maxSocketNameLength = 0;
  nodeData.inputs.forEach((input: NodeInput) => {
    let socket = new ClassicPreset.Socket(input.name);
    if (!node.inputs.hasOwnProperty(input.name)) {
      const inp = new ClassicPreset.Input(socket, input.name);
      node.addInput(input.name, inp);
      maxSocketNameLength = Math.max(maxSocketNameLength, input.name.length);
    }
  });

  nodeData.outputs.forEach((output: NodeOutput) => {
    let socket = new ClassicPreset.Socket(output.name);
    if (!node.outputs.hasOwnProperty(output.name)) {
      node.addOutput(output.name, new ClassicPreset.Output(socket, output.name));
      maxSocketNameLength = Math.max(maxSocketNameLength, output.name.length);
    }
  });
  node.height = Math.max(140, node.height + (nodeData.inputs.length + nodeData.outputs.length) * 35)
  node.width += maxSocketNameLength * 5;

  return node;
}

const customPadding = () => ({
  top: 80,
  left: 30,
  right: 30,
  bottom: 50
});

export async function createEditor(container: HTMLElement, workgraphData: any) {
  container.innerHTML = ''

  const editor = new NodeEditor<Schemes>();
  const area = new AreaPlugin<Schemes, AreaExtra>(container);
  const connection = new ConnectionPlugin<Schemes, AreaExtra>();
  const render = new ReactPlugin<Schemes, AreaExtra>();
  const scopes = new ScopesPlugin<Schemes>({padding: customPadding});
  const arrange = new AutoArrangePlugin<Schemes>();
  const contextMenu = new ContextMenuPlugin<Schemes>({
    items: ContextMenuPresets.classic.setup([
    ])
  });
  const minimap = new MinimapPlugin<Schemes>({
    boundViewport: true
  });

  AreaExtensions.selectableNodes(area, AreaExtensions.selector(), {
    accumulating: AreaExtensions.accumulateOnCtrl()
  });
  AreaExtensions.showInputControl(area);

  render.addPreset(Presets.classic.setup());
  render.addPreset(Presets.contextMenu.setup());
  render.addPreset(Presets.minimap.setup({ size: 200 }));
  render.addPreset(
    Presets.classic.setup({
      customize: {
        control(data) {
          if (data.payload instanceof AtomsControl) {
            return AtomsItem;
          }
          if (data.payload instanceof ClassicPreset.InputControl) {
            return Presets.classic.Control;
          }
          return null;
        }
      }
    })
  );

  connection.addPreset(ConnectionPresets.classic.setup());
  scopes.addPreset(ScopesPresets.classic.setup());

  const applier = new ArrangeAppliers.TransitionApplier<Schemes, never>({
    duration: 500,
    timingFunction: (t) => t,
    async onTick() {
      await AreaExtensions.zoomAt(area, editor.getNodes());
    }
  });

  arrange.addPreset(ArrangePresets.classic.setup());

  editor.use(area);
  // area.use(connection);
  area.use(render);
  area.use(scopes);
  area.use(arrange);
  area.use(contextMenu);
  area.use(minimap);

  await loadJSON(editor, area, workgraphData);

  async function layout(animate: boolean) {
    await arrange.layout({ applier: animate ? applier : undefined });
    AreaExtensions.zoomAt(area, editor.getNodes());
  };

  await layout(true)

  return {
    editor: editor,
    area: area,
    layout: layout,
    destroy: () => area.destroy()
  };
}
