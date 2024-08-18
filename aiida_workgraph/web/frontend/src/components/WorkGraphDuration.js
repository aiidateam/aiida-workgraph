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
    const [useItemType, setUseItemType] = useState("task");

    const fetchData = async () => {
        try {
            const response = await fetch(`http://localhost:8000/api/workgraph-state/${id}?item_type=${useItemType}`);
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
        setInitialLoad(true);
        fetchData();
        const interval = setInterval(fetchData, 5000);
        return () => clearInterval(interval);
    }, [id, useItemType]);

    useEffect(() => {
        if (Object.keys(processesInfo).length) {
            const newGroups = Object.keys(processesInfo).map((key, idx) => ({
                id: idx,
                title: key
            }));

            const newItems = Object.entries(processesInfo).map(([key, { ctime, mtime, process_type }], idx) => ({
                id: idx,
                group: idx,
                title: key,
                start_time: ctime ? moment(ctime) : null,
                end_time: mtime ? moment(mtime) : null
            }));

            setGroups(newGroups);
            setItems(newItems);

            if (initialLoad) {
                const validStartTimes = newItems.map(item => item.start_time).filter(time => time);
                const validEndTimes = newItems.map(item => item.end_time).filter(time => time);
                if (validStartTimes.length && validEndTimes.length) {
                    setTimeStart(moment.min(validStartTimes).valueOf());
                    setTimeEnd(moment.max(validEndTimes).valueOf());
                }
                else {
                    // use the current time as the start time
                    setTimeStart(moment().valueOf());
                    // use the current time + 1 hour as the end time
                    setTimeEnd(moment().add(1, 'hour').valueOf());
                }
                setInitialLoad(false);
            }
        } else {
            setGroups([]);
            setItems([]);
            setTimeStart(null);
            setTimeEnd(null);
        }
    }, [processesInfo]);

    const minZoom = 10000; // 10 seconds in milliseconds
    const maxZoom = 365.25 * 24 * 60 * 60 * 1000; // 1 year in milliseconds

    return (
        <div style={{ padding: '10px', margin: '20px', border: '1px solid #ccc' }}>
            <h1 style={{ textAlign: 'center', color: '#2a3f5f' }}>Node Process Timeline</h1>
            <div style={{textAlign: 'left', fontSize: '16px', color: '#555' }}>
                <p>
                The timeline uses bars to represent the active periods of process nodes, marked from their creation (ctime) to their last modification (mtime). It's important to note that these timestamps do not necessarily correlate with the actual running time since processes might be queued or paused.
                </p>
                <p>
                    Data nodes are shown as rows without bars.
                </p>
            </div>
            <div style={{ margin: '10px' }}>
    <label>
        Task:
        <input
            type="radio"
            value="task"
            checked={useItemType === "task"}
            onChange={(e) => setUseItemType(e.target.value)}
        />
    </label>
    <label style={{ marginLeft: '20px' }}>
        Called Process:
        <input
            type="radio"
            value="called_process"
            checked={useItemType === "called_process"}
            onChange={(e) => setUseItemType(e.target.value)}
        />
    </label>
</div>
            {items.length > 0 ? (
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
            ) : (
                <div style={{ textAlign: 'center', color: '#D32F2F', marginTop: '20px', fontSize: '18px', fontWeight: 'bold' }}>
                    There are no items to display.
                </div>
            )}
        </div>
    );
};

export default NodeDurationGraph;
