import React, { useState, useEffect } from 'react';
import ReactPaginate from 'react-paginate';
import { Link } from 'react-router-dom'; // Import Link


function WorkTree() {
    const [data, setData] = useState([]);
    const [currentPage, setCurrentPage] = useState(0); // `react-paginate` uses zero-based indexing
    const [itemsPerPage] = useState(10);

    useEffect(() => {
        // Replace with your actual endpoint
        fetch('http://localhost:8000/worktree-data')
            .then(response => response.json())
            .then(data => setData(data))
            .catch(error => console.error('Error fetching data: ', error));
    }, []);

    const indexOfLastItem = (currentPage + 1) * itemsPerPage;
    const indexOfFirstItem = indexOfLastItem - itemsPerPage;
    const currentItems = data.slice(indexOfFirstItem, indexOfLastItem);

    const handlePageClick = (event) => {
        setCurrentPage(event.selected);
    };

    return (
        <div>
            <h2>WorkTree</h2>
            <table className="table">
                <thead>
                    <tr>
                        <th>PK</th>
                        <th>Created</th>
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
                // Add these if needed for additional styling
                previousClassName={'previousButton'}
                nextClassName={'nextButton'}
                breakClassName={'pageBreak'}
            />
        </div>
    );
}

export default WorkTree;
