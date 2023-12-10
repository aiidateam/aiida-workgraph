import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';

function WorkTreeDetail() {
    let { pk } = useParams();
    const [worktreeData, setWorktreeData] = useState(null);

    useEffect(() => {
        fetch(`http://localhost:8000/worktree/${pk}`)
            .then(response => response.json())
            .then(data => setWorktreeData(data))
            .catch(error => console.error('Error fetching data:', error));
    }, [pk]); // useEffect will re-run if pk changes

    return (
        <div>
            <h2>WorkTree Detail - PK: {pk}</h2>
            {/* Render worktreeData here */}
            {worktreeData && (
                <div>
                    {/* Display the data. For example: */}
                    <p>State: {worktreeData.state}</p>
                    <p>Create Time: {worktreeData.ctime}</p>
                    {/* Render more data as needed */}
                </div>
            )}
        </div>
    );
}

export default WorkTreeDetail;
