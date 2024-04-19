import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';

function WorkGraphDetail() {
    let { pk } = useParams();
    const [workgraphData, setWorktreeData] = useState(null);

    useEffect(() => {
        fetch(`http://localhost:8000/workgraph/${pk}`)
            .then(response => response.json())
            .then(data => setWorktreeData(data))
            .catch(error => console.error('Error fetching data:', error));
    }, [pk]); // useEffect will re-run if pk changes

    return (
        <div>
            <h2>WorkGraph Detail - PK: {pk}</h2>
            {/* Render workgraphData here */}
            {workgraphData && (
                <div>
                    {/* Display the data. For example: */}
                    <p>State: {workgraphData.state}</p>
                    <p>Create Time: {workgraphData.ctime}</p>
                    {/* Render more data as needed */}
                </div>
            )}
        </div>
    );
}

export default WorkGraphDetail;
