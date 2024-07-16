html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Rete.js with React in Vanilla JS</title>
    <!-- Import React, ReactDOM, and Babel from CDN -->
    <script src="https://unpkg.com/react@18.2.0/umd/react.development.js" crossorigin></script>
    <script src="https://unpkg.com/react-dom@18.2.0/umd/react-dom.development.js" crossorigin></script>
    <script src="https://unpkg.com/@babel/standalone/babel.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/react-is/18.2.0/umd/react-is.production.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/styled-components@5.3.6/dist/styled-components.min.js"></script>
    <script src="https://unpkg.com/elkjs@0.8.2/lib/elk.bundled.js"></script>

    <!-- Import Rete.js and its plugins from CDN -->
    <script src="https://cdn.jsdelivr.net/npm/rete@2.0.3/rete.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-area-plugin@2.0.3/rete-area-plugin.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-connection-plugin@2.0.2/rete-connection-plugin.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-render-utils@2.0.2/rete-render-utils.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-react-plugin@2.0.5/rete-react-plugin.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-auto-arrange-plugin@2.0.1/rete-auto-arrange-plugin.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/rete-minimap-plugin@2.0.1/rete-minimap-plugin.min.js"></script>

    <style>
        .App {
            font-family: sans-serif;
            background: rgb(200, 190, 190);
        }
        .rete {
          position: relative;
          font-size: 1rem;
          margin: 1em;
          border-radius: 1em;
          text-align: left;
        }
        #fullscreen-btn {
            margin-left: 10px;
        }
        body {
            overflow: hidden;
            margin: 0;
            padding: 0;
        }
    </style>
</head>
<body>
    <div id="root"></div>
    <script type="text/babel">

        const { useState, useRef, useEffect } = React;
        const { createRoot } = ReactDOM;
        const { NodeEditor, ClassicPreset } = Rete;
        const { AreaPlugin, AreaExtensions } = ReteAreaPlugin;
        const { ConnectionPlugin, Presets: ConnectionPresets } = ReteConnectionPlugin;
        const { ReactPlugin, Presets } = ReteReactPlugin;
        const { AutoArrangePlugin, Presets: ArrangePresets, ArrangeAppliers} = ReteAutoArrangePlugin;
        const { MinimapExtra, MinimapPlugin } = ReteMinimapPlugin;
        const { RenderUtils } = ReteRenderUtils;
        const styled = window.styled;

        const workgraphData = __WORKGRAPH_DATA__

        // Define Schemes to use in vanilla JS
        const Schemes = {
            Node: ClassicPreset.Node,
            Connection: ClassicPreset.Connection
        };

        class Node extends ClassicPreset.Node {
          width = 180;
          height = 100;
        }
        class Connection extends ClassicPreset.Connection {}

        function createDynamicNode(nodeData) {
          const node = new Node(nodeData.label);
          // resize the node based on the max length of the input/output names
          let maxSocketNameLength = 0;
          nodeData.inputs.forEach((input) => {
            let socket = new ClassicPreset.Socket(input.name);
            if (!node.inputs.hasOwnProperty(input.name)) {
              node.addInput(input.name, new ClassicPreset.Input(socket, input.name));
              maxSocketNameLength = Math.max(maxSocketNameLength, input.name.length);
            }
          });

          nodeData.outputs.forEach((output) => {
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


        async function addNode(editor, area, nodeData) {
          console.log("Adding node", nodeData);
          const node = createDynamicNode(nodeData);
          await editor.addNode(node);
          editor.nodeMap[nodeData.label] = node; // Assuming each nodeData has a unique ID
          await area.translate(node.id, { x: nodeData.position[0], y: nodeData.position[1] });
        }

        async function addLink(editor, area, layout, linkData) {
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

        async function loadJSON(editor, area, layout, workgraphData) {
          for (const nodeId in workgraphData.nodes) {
            const nodeData = workgraphData.nodes[nodeId];
            await addNode(editor, area, nodeData);
          }

          // Adding connections based on workgraphData
          workgraphData.links.forEach(async (link) => { // Specify the type of link here
            await addLink(editor, area, layout, link);
          });
        }

        async function createEditor(container) {
            const socket = new ClassicPreset.Socket("socket");

            const editor = new NodeEditor(Schemes);
            const area = new AreaPlugin(container);
            const connection = new ConnectionPlugin();
            const render = new ReactPlugin({ createRoot });
            const arrange = new AutoArrangePlugin();

            const minimap = new MinimapPlugin({
              boundViewport: true
            });

            AreaExtensions.selectableNodes(area, AreaExtensions.selector(), {
                accumulating: AreaExtensions.accumulateOnCtrl(),
            });

            render.addPreset(Presets.classic.setup());
            render.addPreset(Presets.minimap.setup({ size: 200 }));

            connection.addPreset(ConnectionPresets.classic.setup());

            const applier = new ArrangeAppliers.TransitionApplier({
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
            area.use(minimap);


            AreaExtensions.simpleNodesOrder(area);

            async function layout(animate) {
              await arrange.layout({ applier: animate ? applier : undefined });
              AreaExtensions.zoomAt(area, editor.getNodes());
            }

            // Adding nodes based on workgraphData
            const nodeMap = {}; // To keep track of created nodes for linking
            editor.nodeMap = nodeMap;


            return {
              editor: editor,
              area: area,
              layout: layout,
              destroy: () => area.destroy()
            };
        }

        function toggleFullScreen() {
            if (!document.fullscreenElement) {
                document.documentElement.requestFullscreen();
            } else if (document.exitFullscreen) {
                document.exitFullscreen();
            }
        }

        function App() {
            const [editor, setEditor] = useState(null);
            const containerRef = useRef(null);

            useEffect(() => {
                if (containerRef.current && !editor) {
                    createEditor(containerRef.current).then((editor) => {
                        setEditor(editor);
                        loadJSON(editor.editor, editor.area, editor.layout, workgraphData).then(() => {
                          // aplly layout twice to ensure all nodes are arranged
                          editor?.layout(false).then(() => editor?.layout(true));
                        });
                    });
                }
                if (document.getElementById('fullscreen-btn')) {
                    document.getElementById('fullscreen-btn').addEventListener('click', toggleFullScreen);
                }
                return () => {
                    if (editor) {
                        editor.destroy();
                    }
                };
            }, [containerRef, editor]);

            return (

                <div className="App">
                    <div>
                      <button onClick={() => editor?.layout(true)}>Arrange</button>
                      <button id="fullscreen-btn">Fullscreen</button>
                    </div>
                    <div ref={containerRef} className="rete" style={{ height: "100vh", width: "100%" }}></div>
                </div>
            );
        }

        const rootElement = document.getElementById("root");
        const root = createRoot(rootElement);

        root.render(
                <App />
        );
    </script>
</body>
</html>
"""
