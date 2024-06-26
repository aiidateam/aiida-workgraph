import React, { useEffect, useState } from 'react';
import Timeline from 'react-calendar-timeline';
import 'react-calendar-timeline/lib/Timeline.css';
import moment from 'moment';

const NodeDurationGraph = ({ id }) => {
    const [processesInfo, setProcessesInfo] = useState({});
    const [groups, setGroups] = useState([]);
    const [items, setItems] = useState([]);

    // Function to fetch data from the backend
    const fetchData = async () => {
        try {
            const response = await fetch(`http://localhost:8000/api/workgraph-state/${id}`);
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            const data = await response.json();
            setProcessesInfo(data);
        } catch (error) {
            console.error('Error fetching data:', error);
        }
    };

    useEffect(() => {
        fetchData();

        // Set up polling (optional)
        const interval = setInterval(fetchData, 1000); // Adjust interval as needed

        return () => clearInterval(interval); // Clear interval on cleanup
    }, [id]);

    useEffect(() => {
        // Transform the processesInfo into groups and items for the timeline
        const newGroups = Object.keys(processesInfo).map((key, idx) => ({
            id: idx,
            title: key
        }));

        const newItems = Object.entries(processesInfo).map(([key, { ctime, mtime }], idx) => ({
            id: idx,
            group: idx,
            title: key,
            start_time: moment(ctime),
            end_time: moment(mtime)
        }));

        setGroups(newGroups);
        setItems(newItems);
    }, [processesInfo]);

    const defaultTimeStart = moment.min(items.map(item => item.start_time));
    const defaultTimeEnd = moment.max(items.map(item => item.end_time));

    // Define minimum and maximum zoom (in milliseconds)
    const minZoom = 10000; // 1 second in milliseconds
    const maxZoom = 365.25 * 24 * 60 * 60 * 1000; // 1 year in milliseconds

    return (
        <Timeline
            groups={groups}
            items={items}
            defaultTimeStart={defaultTimeStart}
            defaultTimeEnd={defaultTimeEnd}
            lineHeight={50}
            minZoom={minZoom}
            maxZoom={maxZoom}
            canMove={false} // Set to false to prevent moving items
            canChangeGroup={false} // Set to false to prevent changing groups
            canResize={'both'} // Allow resizing items from both sides
        />
    );
};

export default NodeDurationGraph;
