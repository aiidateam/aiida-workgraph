import * as React from "react";
import { createRender, useModelState } from "@anywidget/react";
// import { createEditor } from "./editor";
import { createEditor } from "./default_rete";
import "./widget.css";



/* Modify the useRete function to support passing worktree data to createEditor */
export function useRete<T extends { destroy(): void }>(
	create: (el: HTMLElement, data: any) => Promise<T>,
	worktreeData: any
  ) {
	const [container, setContainer] = React.useState<null | HTMLElement>(null);
	const editorRef = React.useRef<T>();
	const [editor, setEditor] = React.useState<T | null>(null);
	const ref = React.useRef(null);

	React.useEffect(() => {
	  if (container) {
		if (editorRef.current) {
		  editorRef.current.destroy();
		  container.innerHTML = '';
		}
		create(container, worktreeData).then((value) => {
		  editorRef.current = value;
		  setEditor(value);
		});
	  }
	}, [container, create, worktreeData]); // Add worktreeData as a dependency

	React.useEffect(() => {
	  return () => {
		if (editorRef.current) {
		  editorRef.current.destroy();
		}
	  };
	}, []);

	React.useEffect(() => {
	  if (ref.current) {
		setContainer(ref.current);
	  }
	}, [ref.current]);

	return [ref, editor] as const;
  }


const render = createRender(() => {
	const [value, setValue] = useModelState<number>("value");
	const [ref] = useRete(createEditor, value);
	return (
		<div className="App">
		  <div ref={ref} className="rete"></div>
		</div>
	  );
});

export default { render };
