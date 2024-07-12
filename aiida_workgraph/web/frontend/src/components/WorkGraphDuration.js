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
        const interval = setInterval(fetchData, 5000); // Adjust interval as needed

        return () => clearInterval(interval); // Clear interval on cleanup
    }, [id]);

    useEffect(() => {
        // Transform the processesInfo into groups and items for the timeline
        const newGroups = Object.keys(processesInfo).map((key, idx) => ({
            id: idx,
            title: key
        }));


        const newItems = Object.entries(processesInfo).map(([key, { ctime, mtime, process_type }], idx) => {
            if (process_type) {
            return {
                id: idx,
                group: idx,
                title: key,
                start_time: ctime ? moment(ctime) : null,
                end_time: mtime ? moment(mtime) : null
            };
            }
            return null;
        }).filter(item => item !== null);

        setGroups(newGroups);
        setItems(newItems);
    }, [processesInfo]);

    const validStartTimes = items.map(item => item.start_time).filter(time => time !== null);
    const validEndTimes = items.map(item => item.end_time).filter(time => time !== null);
    const defaultTimeStart = moment.min(validStartTimes);
    const defaultTimeEnd = moment.max(validEndTimes);

    // Define minimum and maximum zoom (in milliseconds)
    const minZoom = 10000; // 1 second in milliseconds
    const maxZoom = 365.25 * 24 * 60 * 60 * 1000; // 1 year in milliseconds

    return (
        <Timeline
            groups={groups}
            items={items}
            visibleTimeStart={defaultTimeStart}
            visibleTimeEnd={defaultTimeEnd}
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
