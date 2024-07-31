// PromptModal.js
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';


function WorkGraphDeleteNodePrompt(props) {
  return (
      <Modal 
        style={{ display: 'block', position: 'initial' }}
        {...props}
      >
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          Modal heading
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <h4>Centered Modal</h4>
        <p>
          Cras mattis consectetur purus sit amet fermentum. Cras justo odio,
          dapibus ac facilisis in, egestas eget quam. Morbi leo risus, porta ac
          consectetur ac, vestibulum at eros.
        </p>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={props.onHide}>Close</Button>
      </Modal.Footer>
      </Modal>
  );
}
//function WorkGraphDeleteNodePrompt({ openOnInit }) {
//    const [open, setOpen] = useState(false);
//    setOpen(openOnInit)
//
//    const handleYesClick = () => {
//      setOpen(false)
//    };
//    const handleNoClick = () => {
//      setOpen(true);
//    };
//
//    return ( 
//      <div
//        className="modal show"
//        style={{ display: 'block', position: 'initial' }}
//      >
//        <Modal isOpen={open}>
//            <h2>Are you sure you want to delete node ...</h2>
//            <button onClick={handleYesClick}>Yes</button>
//            <button onClick={handleNoClick}>No</button>
//        </Modal>
//      </div>);
//}

export default WorkGraphDeleteNodePrompt;
