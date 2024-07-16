import React, { useEffect, useState } from 'react';
import Timeline from 'react-calendar-timeline';
import 'react-calendar-timeline/lib/Timeline.css';
import moment from 'moment';

const NodeDurationGraph = ({ id }) => {
    const [processesInfo, setProcessesInfo] = useState({});
    const [groups, setGroups] = useState([]);
    const [items, setItems] = useState([]);
    const [timeStart, setTimeStart] = useState(null);
    const [timeEnd, setTimeEnd] = useState(null);
    const [initialLoad, setInitialLoad] = useState(true);

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
        const interval = setInterval(fetchData, 5000);
        return () => clearInterval(interval);
    }, [id]);

    useEffect(() => {
        if (Object.keys(processesInfo).length) {
            const newGroups = Object.keys(processesInfo).map((key, idx) => ({
                id: idx,
                title: key
            }));

            const newItems = Object.entries(processesInfo).map(([key, { ctime, mtime, process_type }], idx) => {
                return process_type ? {
                    id: idx,
                    group: idx,
                    title: key,
                    start_time: ctime ? moment(ctime) : null,
                    end_time: mtime ? moment(mtime) : null
                } : null;
            }).filter(item => item !== null);

            setGroups(newGroups);
            setItems(newItems);

            if (initialLoad) {
                const validStartTimes = newItems.map(item => item.start_time).filter(time => time);
                const validEndTimes = newItems.map(item => item.end_time).filter(time => time);
                if (validStartTimes.length && validEndTimes.length) {
                    setTimeStart(moment.min(validStartTimes).valueOf());
                    setTimeEnd(moment.max(validEndTimes).valueOf());
                }
                setInitialLoad(false);
            }
        }
    }, [processesInfo]);

    const minZoom = 10000; // 10 seconds in milliseconds
    const maxZoom = 365.25 * 24 * 60 * 60 * 1000; // 1 year in milliseconds

    if (!timeStart || !timeEnd) {
        return <div>Loading timeline...</div>;
    }

    return (
            <div style={{ padding: '10px', margin: '20px', border: '1px solid #ccc',  }}>
            <h1 style={{ textAlign: 'center', color: '#2a3f5f' }}>Node Process Timeline</h1>
            <div style={{textAlign: 'left', fontSize: '16px', color: '#555' }}>
                <p>
                The timeline uses bars to represent the active periods of process nodes, marked from their creation (ctime) to their last modification (mtime). It's important to note that these timestamps do not necessarily correlate with the actual running time since processes might be queued or paused.
                </p>
                <p>
                    Data nodes are shown as rows without bars.
                </p>
                </div>
            <Timeline
                groups={groups}
                items={items}
                visibleTimeStart={timeStart}
                visibleTimeEnd={timeEnd}
                onTimeChange={(visibleTimeStart, visibleTimeEnd, updateScrollCanvas) => {
                    setTimeStart(visibleTimeStart);
                    setTimeEnd(visibleTimeEnd);
                    updateScrollCanvas(visibleTimeStart, visibleTimeEnd);
                }}
                lineHeight={50}
                minZoom={minZoom}
                maxZoom={maxZoom}
                canMove={false}
                canChangeGroup={false}
                canResize={'both'}
            />
        </div>
    );
};

export default NodeDurationGraph;
