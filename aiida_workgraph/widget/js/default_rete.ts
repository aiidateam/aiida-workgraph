import { NodeEditor, GetSchemes, ClassicPreset } from "rete";
import { AreaPlugin, AreaExtensions } from "rete-area-plugin";
import {
  ConnectionPlugin,
  Presets as ConnectionPresets
} from "rete-connection-plugin";
import { ReactPlugin, Presets, ReactArea2D } from "rete-react-plugin";
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

// May move these interfaces to a separate file
interface NodeInput {
  name: string;
}

interface NodeOutput {
  name: string;
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
}

interface NodeMap {
  [key: string]: Node;
}



class Node extends ClassicPreset.Node {
  width = 180;
  height = 100;
}
class Connection<N extends Node> extends ClassicPreset.Connection<N, N> {}

type Schemes = GetSchemes<Node, Connection<Node>>;
type AreaExtra = ReactArea2D<any> | MinimapExtra | ContextMenuExtra;


export function createDynamicNode(nodeData: any) {
  const node = new Node(nodeData.label);
  // resize the node based on the max length of the input/output names
  let maxSocketNameLength = 0;
  nodeData.inputs.forEach((input: NodeInput) => {
    let socket = new ClassicPreset.Socket(input.name);
    if (!node.inputs.hasOwnProperty(input.name)) {
      node.addInput(input.name, new ClassicPreset.Input(socket, input.name));
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

export async function loadJSON(editor, area, layout, workgraphData) {
  for (const nodeId in workgraphData.nodes) {
    const nodeData = workgraphData.nodes[nodeId];
    await addNode(editor, area, nodeData);
  }

  // Adding connections based on workgraphData
  workgraphData.links.forEach(async (link: LinkData) => { // Specify the type of link here
    await addLink(editor, area, layout, link);
  });
}

export async function addNode(editor, area, nodeData) {
  console.log("Adding node", nodeData);
  const node = createDynamicNode(nodeData);
  await editor.addNode(node);
  editor.nodeMap[nodeData.label] = node; // Assuming each nodeData has a unique ID
  await area.translate(node.id, { x: nodeData.position[0], y: nodeData.position[1] });
}

export async function addLink(editor, area, layout, linkData) {
  const fromNode = editor.nodeMap[linkData.from_node];
  const toNode = editor.nodeMap[linkData.to_node];
  console.log("fromNode", fromNode, "toNode", toNode);
  let socket;
  if (fromNode && toNode) {
    socket = new ClassicPreset.Socket(linkData.from_socket);
    if (!fromNode.outputs.hasOwnProperty(linkData.from_socket)) {
      fromNode.addOutput(linkData.from_socket, new ClassicPreset.Output(socket, linkData.from_socket));
      fromNode.height += 25; // Increase height of node for each output
      area.update('node', fromNode.id);
    }
    socket = new ClassicPreset.Socket(linkData.to_socket);
    if (!toNode.inputs.hasOwnProperty(linkData.to_socket)) {
      toNode.addInput(linkData.to_socket, new ClassicPreset.Input(socket, linkData.to_socket));
      toNode.height += 25; // Increase height of node for each input
      area.update('node', toNode.id);
    }
    await editor.addConnection(new Connection(fromNode, linkData.from_socket, toNode, linkData.to_socket));
    // await layout(true);

  }
}

export async function removeLink(editor, linkData) {
  const connections = editor.getConnections();
  const connection = connections.find((connection) => {
    return connection.source === editor.nodeMap[linkData.from_node].id && connection.target === editor.nodeMap[linkData.to_node].id;
  });
  if (connection) {
    editor.removeConnection(connection.id);
    console.log("Deleted link successfully");
  }
}

export async function removeNode(editor, name) {
  // remove all connections to the node
  const node = editor.nodeMap[name];
  const connections = editor.getConnections()

  const incomingConnections = connections.filter(connection => connection.target === node.id)
  const outgoingConnections = connections.filter(connection => connection.source === node.id)
  incomingConnections.forEach(connection => {
    editor.removeConnection(connection.id)
  })
  outgoingConnections.forEach(connection => {
    editor.removeConnection(connection.id)
  })
  editor.removeNode(editor.nodeMap[name].id).then(() => {
    delete editor.nodeMap[name];
    console.log("Deleted node successfully");
  });
}

export async function createEditor(container: HTMLElement, settings: any) {
  container.innerHTML = ''

  const editor = new NodeEditor<Schemes>();
  const area = new AreaPlugin<Schemes, AreaExtra>(container);
  const connection = new ConnectionPlugin<Schemes, AreaExtra>();
  const render = new ReactPlugin<Schemes, AreaExtra>();
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

  render.addPreset(Presets.classic.setup());
  render.addPreset(Presets.contextMenu.setup());
  render.addPreset(Presets.minimap.setup({ size: 200 }));

  connection.addPreset(ConnectionPresets.classic.setup());

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
  area.use(arrange);
  area.use(contextMenu);
  if (settings.minimap) {
    area.use(minimap);
  }

  AreaExtensions.simpleNodesOrder(area);


  async function layout(animate: boolean) {
    await arrange.layout({ applier: animate ? applier : undefined });
    AreaExtensions.zoomAt(area, editor.getNodes());
  };

  // Adding nodes based on workgraphData
  const nodeMap: NodeMap = {}; // To keep track of created nodes for linking
  editor.nodeMap = nodeMap;

  // aplly layout twice to ensure all nodes are arranged
  // await layout(true);
  // await layout(true);

  return {
    editor: editor,
    area: area,
    layout: layout,
    destroy: () => area.destroy()
  };
}
