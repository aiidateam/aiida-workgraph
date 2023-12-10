import React, { useState, useEffect } from 'react';
import ReactPaginate from 'react-paginate';
import { Link } from 'react-router-dom';

function WorkTree() {
    const [data, setData] = useState([]);
    const [currentPage, setCurrentPage] = useState(0);
    const [itemsPerPage] = useState(15);
    const [sortField, setSortField] = useState("pk");
    const [sortOrder, setSortOrder] = useState('asc'); // 'asc' or 'desc'
    const [searchQuery, setSearchQuery] = useState('');

    useEffect(() => {
        fetch(`http://localhost:8000/worktree-data?search=${searchQuery}`)
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

    return (
        <div>
            <h2>WorkTree</h2>
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
                        <th onClick={() => sortData('ctime')}>Created</th>
                        <th>Process Label</th>
                        <th>State</th>
                    </tr>
                </thead>
                <tbody>
                    {currentItems.map(item => (
                        <tr key={item.pk}>
                            <td>
                                <Link to={`/worktree/${item.pk}`}>{item.pk}</Link>
                            </td>
                            <td>{item.ctime}</td>
                            <td>{item.process_label}</td>
                            <td>{item.state}</td>
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
                // Additional styling classes
                previousClassName={'previousButton'}
                nextClassName={'nextButton'}
                breakClassName={'pageBreak'}
            />
        </div>
    );
}

export default WorkTree;
