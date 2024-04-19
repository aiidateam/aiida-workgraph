import React, { useState, useEffect } from 'react';
import ReactPaginate from 'react-paginate';
import { Link } from 'react-router-dom';
import { toast, ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import { FaPlay, FaPause, FaTrash } from 'react-icons/fa';
import './WorkGraphTable.css';

function DataNode() {
    const [data, setData] = useState([]);
    const [currentPage, setCurrentPage] = useState(0);
    const [itemsPerPage] = useState(15);
    const [sortField, setSortField] = useState("pk");
    const [sortOrder, setSortOrder] = useState('desc');
    const [searchTypeQuery, setSearchTypeQuery] = useState('');
    const [searchLabelQuery, setSearchLabelQuery] = useState(''); // New state for Label search

    useEffect(() => {
        fetch(`http://localhost:8000/api/datanode-data?typeSearch=${searchTypeQuery}&labelSearch=${searchLabelQuery}`) // Include labelSearch query parameter
            .then(response => response.json())
            .then(data => setData(data))
            .catch(error => console.error('Error fetching data: ', error));
    }, [searchTypeQuery, searchLabelQuery]); // Include searchLabelQuery as a dependency

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

    const handleDeleteClick = (item) => {
        fetch(`http://localhost:8000/api/workgraph/delete/${item.pk}`, {
            method: 'DELETE',
        })
        .then(response => response.json())
        .then(data => {
            if (data.message) {
                toast.success(data.message);
                fetch(`http://localhost:8000/api/workgraph-data?TypeSearch=${searchTypeQuery}&labelSearch=${searchLabelQuery}`) // Include labelSearch in refresh
                    .then(response => response.json())
                    .then(data => setData(data))
                    .catch(error => console.error('Error fetching data: ', error));
            } else {
                toast.error('Error deleting item');
            }
        })
        .catch(error => console.error('Error deleting item: ', error));
    };

    return (
        <div>
            <h2>Data Node</h2>
            <div className="search-container">
                <input
                    type="text"
                    placeholder="Search by Type..."
                    value={searchTypeQuery}
                    onChange={(e) => setSearchTypeQuery(e.target.value)}
                    className="search-input"
                />
                <input
                    type="text"
                    placeholder="Search by Label..." // Add placeholder for Label
                    value={searchLabelQuery} // Create state for Label search
                    onChange={(e) => setSearchLabelQuery(e.target.value)} // Create handler for Label search
                    className="search-input"
                />
            </div>
            <table className="table">
                <thead>
                    <tr>
                        <th onClick={() => sortData('pk')}>PK {renderSortIndicator('pk')}</th>
                        <th onClick={() => sortData('ctime')}>Created {renderSortIndicator('ctime')}</th>
                        <th>Type</th>
                        <th>Label</th>
                        <th>Actions</th>
                    </tr>
                </thead>
                <tbody>
                    {currentItems.map(item => (
                        <tr key={item.pk}>
                            <td>
                                <Link to={`/datanode/${item.pk}`}>{item.pk}</Link>
                            </td>
                            <td>{item.ctime}</td>
                            <td>{item.node_type}</td>
                            <td>{item.label}</td>
                            <td>
                                <button onClick={() => handleDeleteClick(item)} className="action-button delete-button"><FaTrash /></button>
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
        </div>
    );
}

export default DataNode;
