import React, { useState, useEffect } from 'react';
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

function Settings() {
  const [workers, setWorkers] = useState([]);

  const fetchWorkers = () => {
    fetch('http://localhost:8000/daemon/worker')
      .then(response => response.json())
      .then(data => setWorkers(Object.values(data)))
      .catch(error => console.error('Failed to fetch workers:', error));
  };

  useEffect(() => {
    fetchWorkers();
    const interval = setInterval(fetchWorkers, 1000); // Poll every 5 seconds
    return () => clearInterval(interval); // Clear interval on component unmount
  }, []);

  const handleDaemonControl = (action) => {
    fetch(`http://localhost:8000/daemon/${action}`, { method: 'POST' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`Daemon operation failed: ${response.statusText}`);
        }
        return response.json();
      })
      .then(data => {
        toast.success(`Daemon ${action}ed successfully`);
        fetchWorkers();
      })
      .catch(error => toast.error(error.message));
  };

  const adjustWorkers = (action) => {
    fetch(`http://localhost:8000/daemon/${action}`, { method: 'POST' })
      .then(response => {
        if (!response.ok) {
          throw new Error(`Failed to ${action} workers: ${response.statusText}`);
        }
        return response.json();
      })
      .then(data => {
        toast.success(`Workers ${action}ed successfully`);
        fetchWorkers(); // Refetch workers after adjusting
      })
      .catch(error => toast.error(error.message));
  };

  return (
    <div>
      <h2>Daemon Control</h2>
      <ToastContainer />
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
          {workers.map(worker => (
            <tr key={worker.pid}>
              <td>{worker.pid}</td>
              <td>{worker.mem}</td>
              <td>{worker.cpu}</td>
              <td>{new Date(worker.started * 1000).toLocaleString()}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <button className="button button-start" onClick={() => handleDaemonControl('start')}>Start Daemon</button>
      <button className="button button-stop" onClick={() => handleDaemonControl('stop')}>Stop Daemon</button>
      <button className="button button-adjust" onClick={() => adjustWorkers('increase')}>Increase Workers</button>
      <button className="button button-adjust" onClick={() => adjustWorkers('decrease')}>Decrease Workers</button>
    </div>
  );
}

export default Settings;
