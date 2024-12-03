import React, { useState, useEffect } from 'react';
import ReactPaginate from 'react-paginate';
import { Link } from 'react-router-dom';
import { toast, ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import { FaPlay, FaPause, FaTrash } from 'react-icons/fa'; // Import icons from react-icons
import WorkGraphConfirmModal from './WorkGraphModals';
import './WorkGraphTable.css'; // Import a custom CSS file for styling


function WorkGraph() {
    const [data, setData] = useState([]);
    const [currentPage, setCurrentPage] = useState(0);
    const [itemsPerPage] = useState(15);
    const [sortField, setSortField] = useState("pk");
    const [sortOrder, setSortOrder] = useState('desc'); // 'asc' or 'desc'
    const [searchQuery, setSearchQuery] = useState('');

    useEffect(() => {
        fetch(`http://localhost:8000/api/workgraph-data?search=${searchQuery}`)
            .then(response => response.json())
            .then(data => setData(data))
            .catch(error => console.error('Error fetching data: ', error));
    }, [searchQuery]);

    const sortData = (field) => {
        const order = field === sortField && sortOrder === 'asc' ? 'desc' : 'asc';
        const sortedData = [...data].sort((a, b) => {
            if (a[field] < b[field]) return order === 'asc' ? -1 : 1;
            if (a[field] > b[field]) return order === 'asc' ? 1 : -1;
            return 0;
        });
        setData(sortedData);
        setSortField(field);
        setSortOrder(order);
    };

    // Function to render sort indicator
    const renderSortIndicator = (field) => {
        if (sortField === field) {
            return sortOrder === 'asc' ? ' 🔼' : ' 🔽';
        }
        return '';
    };

    const indexOfLastItem = (currentPage + 1) * itemsPerPage;
    const indexOfFirstItem = indexOfLastItem - itemsPerPage;
    const currentItems = data.slice(indexOfFirstItem, indexOfLastItem);

    const handlePageClick = (event) => {
        setCurrentPage(event.selected);
    };

    // Function to handle pause click
    const handlePauseClick = (item) => {
        // Make an API request to pause the workgraph item
        fetch(`http://localhost:8000/api/workgraph/pause/${item.pk}`, {
            method: 'POST',
        })
        .then(response => response.json())
        .then(data => {
            // Show toast notification for success or error
            if (data.message) {
                toast.success(data.message);
                // Refresh the table after pause (you may fetch the updated data here)
                fetch(`http://localhost:8000/api/workgraph-data?search=${searchQuery}`)
                    .then(response => response.json())
                    .then(data => setData(data))
                    .catch(error => console.error('Error fetching data: ', error));
            } else {
                toast.error('Error pausing item');
            }
        })
        .catch(error => console.error('Error pausing item: ', error));
    };

    // Function to handle play click
    const handlePlayClick = (item) => {
        // Make an API request to play the workgraph item
        fetch(`http://localhost:8000/api/workgraph/play/${item.pk}`, {
            method: 'POST',
        })
        .then(response => response.json())
        .then(data => {
            // Show toast notification for success or error
            if (data.message) {
                toast.success(data.message);
                // Refresh the table after play (you may fetch the updated data here)
                fetch(`http://localhost:8000/api/workgraph-data?search=${searchQuery}`)
                    .then(response => response.json())
                    .then(data => setData(data))
                    .catch(error => console.error('Error fetching data: ', error));
            } else {
                toast.error('Error playing item');
            }
        })
        .catch(error => console.error('Error playing item: ', error));
    };

    // Function to handle delete click
      const handleDeleteNode = (item) => {
        // Make an API request to delete the workgraph item
        fetch(`http://localhost:8000/api/workgraph/delete/${item.pk}`, {
            method: 'DELETE',
        })
        .then(response => response.json())
        .then(data => {
            // Show toast notification for success or error
            if (data.deleted) {
                toast.success(data.message);
                // Refresh the table after delete (you may fetch the updated data here)
                fetch(`http://localhost:8000/api/workgraph-data?search=${searchQuery}`)
                    .then(response => response.json())
                    .then(data => setData(data))
                    .catch(error => console.error('Error fetching data: ', error));
            } else {
                toast.error('Error deleting item');
            }
        })
        .catch(error => console.error('Error deleting item: ', error));
    };

    const [toDeleteItem, setToDeleteItem] = useState(null);
    const [showConfirmDeleteModal, setShowConfirmDeleteModal] = useState(false);
    const [bodyTextConfirmDeleteModal, setBodyTextConfirmDeleteModal] =  useState(<p></p>);
    // useEffect to ensure this happens after the toDeletItem has been updated by the delete button click
    useEffect(() => {
        if (toDeleteItem != null) {
          fetch(`http://localhost:8000/api/workgraph/delete/${toDeleteItem.pk}?dry_run=True`, {
              method: 'DELETE',
          })
          .then(response => response.json())
          .then(data => {
              if (data.deleted_nodes.length > 0) {
                // the node {item.pk} will be always in the list if the deletion is successfull
                if (data.deleted_nodes.length > 1) {
                  let formatted_pks = data.deleted_nodes.map((x) => ` ${x.toString()}`);
                  data.deleted_nodes.splice(data.deleted_nodes.indexOf(toDeleteItem.pk), 1)
                  setBodyTextConfirmDeleteModal(
                    <p>
                    Are you sure you want to delete node PK&lt;{toDeleteItem.pk}&gt; and {data.deleted_nodes.length} its dependent nodes?
                    <b> A deletion is irreversible.</b>
                    <br/><br/>
                    List of PKs of nodes dependent on the node PK&lt;{toDeleteItem.pk}&gt; that also will be deleted:
                    <br/> {formatted_pks.toString()}
                    </p>
                  );
                } else {
                  setBodyTextConfirmDeleteModal(
                    <p>
                    Are you sure you want to delete node {toDeleteItem.pk}?
                    <b> A deletion is irreversible.</b>
                    </p>
                  );
                }
                setShowConfirmDeleteModal(true)
              } else {
                  toast.error('Error deleting item.');
              }
          })
          .catch(error => console.error('Error deleting item: ', error));
        };
      }, [toDeleteItem]);

    return (
        <div>
            <h2>WorkGraph</h2>
            <div className="search-container">
                <input
                    type="text"
                    placeholder="Search..."
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    className="search-input"
                />
            </div>
            <table className="table">
                <thead>
                    <tr>
                        <th onClick={() => sortData('pk')}>PK {renderSortIndicator('pk')}</th>
                        <th onClick={() => sortData('ctime')}>Created {renderSortIndicator('ctime')}</th>
                        <th>Process Label</th>
                        <th>State</th>
                        <th>Actions</th>
                    </tr>
                </thead>
                <tbody>
                    {currentItems.map(item => (
                        <tr key={item.pk}>
                            <td>
                                <Link to={`/workgraph/${item.pk}`}>{item.pk}</Link>
                            </td>
                            <td>{item.ctime}</td>
                            <td>{item.process_label}</td>
                            <td>{item.state}</td>
                            <td>
                                <button onClick={() => handlePauseClick(item)} className="action-button pause-button"><FaPause /></button>
                                <button onClick={() => handlePlayClick(item)} className="action-button play-button"><FaPlay /></button>
                                <button
                                  onClick={
                                    () =>
                                    {
                                      // we need to copy it to change the reference
                                      // so useEffect hook gets activated
                                      setToDeleteItem(structuredClone(item));
                                    }
                                  }
                                  className="action-button delete-button"><FaTrash />
                                </button>
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>
            <ReactPaginate
                previousLabel={'previous'}
                nextLabel={'next'}
                breakLabel={'...'}
                pageCount={Math.ceil(data.length / itemsPerPage)}
                marginPagesDisplayed={2}
                pageRangeDisplayed={5}
                onPageChange={handlePageClick}
                containerClassName={'pagination'}
                activeClassName={'active'}
                previousClassName={'previousButton'}
                nextClassName={'nextButton'}
                breakClassName={'pageBreak'}
            />
            <ToastContainer autoClose={3000} />
            <WorkGraphConfirmModal
                show={showConfirmDeleteModal}
                setShow={setShowConfirmDeleteModal}
                confirmAction={() => handleDeleteNode(toDeleteItem)}
                cancelAction={() => {}}
                bodyText={bodyTextConfirmDeleteModal}
            />
        </div>
    );
}

export default WorkGraph;
