import React, { useState, useEffect } from 'react';
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

function Settings() {
  const [taskWorkers, setTaskWorkers] = useState([]);
  const [schedulerWorkers, setSchedulerWorkers] = useState([]);

  // Fetching task workers
  const fetchTaskWorkers = () => {
    fetch('http://localhost:8000/api/daemon/task/worker')
      .then(response => response.json())
      .then(data => setTaskWorkers(Object.values(data)))
      .catch(error => console.error('Failed to fetch task workers:', error));
  };

  // Fetching scheduler workers
  const fetchSchedulerWorkers = () => {
    fetch('http://localhost:8000/api/daemon/scheduler/worker')
      .then(response => response.json())
      .then(data => setSchedulerWorkers(Object.values(data)))
      .catch(error => console.error('Failed to fetch scheduler workers:', error));
  };

  useEffect(() => {
    fetchTaskWorkers();
    fetchSchedulerWorkers();
    const taskInterval = setInterval(fetchTaskWorkers, 1000);
    const schedulerInterval = setInterval(fetchSchedulerWorkers, 1000);
    return () => {
      clearInterval(taskInterval);
      clearInterval(schedulerInterval);
    }; // Clear intervals on component unmount
  }, []);

  const handleDaemonControl = (daemonType, action) => {
    fetch(`http://localhost:8000/api/daemon/${daemonType}/${action}`, { method: 'POST' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`${daemonType} daemon operation failed: ${response.statusText}`);
        }
        return response.json();
      })
      .then(() => {
        toast.success(`${daemonType} daemon ${action}ed successfully`);
        if (daemonType === 'task') {
          fetchTaskWorkers();
        } else {
          fetchSchedulerWorkers();
        }
      })
      .catch(error => toast.error(error.message));
  };

  const adjustWorkers = (daemonType, action) => {
    fetch(`http://localhost:8000/api/daemon/${daemonType}/${action}`, { method: 'POST' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`Failed to ${action} workers for ${daemonType}: ${response.statusText}`);
        }
        return response.json();
      })
      .then(() => {
        toast.success(`${daemonType} Workers ${action}ed successfully`);
        if (daemonType === 'task') {
          fetchTaskWorkers();
        } else {
          fetchSchedulerWorkers();
        }
      })
      .catch(error => toast.error(error.message));
  };

  return (
    <div>
      <ToastContainer />
      <div>
        <h3>Task Daemon Control</h3>
        <button className="button button-start" onClick={() => handleDaemonControl('task', 'start')}>Start Task Daemon</button>
        <button className="button button-stop" onClick={() => handleDaemonControl('task', 'stop')}>Stop Task Daemon</button>
        <button className="button button-adjust" onClick={() => adjustWorkers('task', 'increase')}>Increase Task Workers</button>
        <button className="button button-adjust" onClick={() => adjustWorkers('task', 'decrease')}>Decrease Task Workers</button>
        <table className="table">
          <thead>
            <tr>
              <th>PID</th>
              <th>Memory %</th>
              <th>CPU %</th>
              <th>Started</th>
            </tr>
          </thead>
          <tbody>
            {taskWorkers.map(worker => (
              <tr key={worker.pid}>
                <td>{worker.pid}</td>
                <td>{worker.mem}</td>
                <td>{worker.cpu}</td>
                <td>{new Date(worker.started * 1000).toLocaleString()}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div>
        <h3>Scheduler Daemon Control</h3>
        <button className="button button-start" onClick={() => handleDaemonControl('scheduler', 'start')}>Start Scheduler Daemon</button>
        <button className="button button-stop" onClick={() => handleDaemonControl('scheduler', 'stop')}>Stop Scheduler Daemon</button>
        <button className="button button-adjust" onClick={() => adjustWorkers('scheduler', 'increase')}>Increase Scheduler Workers</button>
        <button className="button button-adjust" onClick={() => adjustWorkers('scheduler', 'decrease')}>Decrease Scheduler Workers</button>
        <table className="table">
          <thead>
            <tr>
              <th>PID</th>
              <th>Memory %</th>
              <th>CPU %</th>
              <th>Started</th>
            </tr>
          </thead>
          <tbody>
            {schedulerWorkers.map(worker => (
              <tr key={worker.pid}>
                <td>{worker.pid}</td>
                <td>{worker.mem}</td>
                <td>{worker.cpu}</td>
                <td>{new Date(worker.started * 1000).toLocaleString()}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default Settings;
