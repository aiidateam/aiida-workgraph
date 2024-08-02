// PromptModal.js
import 'bootstrap/dist/css/bootstrap.min.css';
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';


function WorkGraphConfirmModal(props) {

  const handleConfirmClick = () => {
    props.confirmAction();
    props.setShow(false);
  };

  const handleCancelClick = () => {
    props.cancelAction();
    props.setShow(false);
  };

  return (
      <Modal 
        show={props.show}
      >
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          Confirm deletion 
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <p>
          {props.bodyText}
          
        </p>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={handleCancelClick}>Cancel</Button>
        <Button onClick={handleConfirmClick}>Confirm</Button>
      </Modal.Footer>
      </Modal>
  );
}
export default WorkGraphDeleteNodePrompt;
