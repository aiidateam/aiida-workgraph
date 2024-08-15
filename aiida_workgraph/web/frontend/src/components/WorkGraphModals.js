// PromptModal.js
import Modal from 'react-bootstrap/Modal';
import Button from 'react-bootstrap/Button';


function WorkGraphConfirmModal({ confirmAction, cancelAction, bodyText, show, setShow }) {

  const handleConfirmClick = () => {
    try {
      confirmAction();
    } catch(error) {
      console.error('While executing the modal confirm action, an action error was raised: ', error)
    } finally {
      setShow(false);
    }
  };

  const handleCancelClick = () => {
    try {
      cancelAction();
    } catch(error) {
      console.error('While executing the modal cancel action, an action error was raised: ', error)
    } finally {
      setShow(false);
    }
  };

  return (
      <Modal
        show={show}
      >
      <Modal.Header closeButton onClick={() => handleCancelClick()}>
        <Modal.Title id="contained-modal-title-vcenter">
          Confirm deletion
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <p>
          {bodyText}
        </p>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={handleCancelClick}>Cancel</Button>
        <Button onClick={handleConfirmClick}>Confirm</Button>
      </Modal.Footer>
      </Modal>
  );
}
export default WorkGraphConfirmModal;
