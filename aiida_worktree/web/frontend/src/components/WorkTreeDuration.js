import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';

const NodeDurationGraph = ({ processesInfo }) => {
    const d3Container = useRef(null);

    useEffect(() => {
        if (processesInfo && d3Container.current) {
            // Clear the existing graph
            d3.select(d3Container.current).selectAll("*").remove();

            // Transform the data
            const data = Object.entries(processesInfo).map(([name, { ctime, mtime }]) => [name, ctime, mtime]);

            const margin = { top: 20, right: 30, bottom: 40, left: 90 };
            const width = 960 - margin.left - margin.right;
            const height = 250 - margin.top - margin.bottom;

            const svg = d3.select(d3Container.current)
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);

            const x = d3.scaleTime()
                .domain([d3.min(data, d => new Date(d[1])), d3.max(data, d => new Date(d[2]))])
                .range([0, width]);

            svg.append("g")
                .attr("transform", `translate(0,${height})`)
                .call(d3.axisBottom(x))
                .selectAll("text")
                .style("font-size", "12px"); // Adjust X-axis label font size here

            const y = d3.scaleBand()
                .range([0, height])
                .domain(data.map(d => d[0]))
                .padding(0.2);

            svg.append("g")
                .call(d3.axisLeft(y))
                .selectAll("text")
                .style("font-size", "12px"); // Adjust Y-axis label font size here

            svg.selectAll("myRect")
                .data(data)
                .join("rect")
                .attr("x", d => x(new Date(d[1])))
                .attr("y", d => y(d[0]))
                .attr("width", d => x(new Date(d[2])) - x(new Date(d[1])))
                .attr("height", y.bandwidth())
                .attr("fill", "#69b3a2");
        }
    }, [processesInfo]);

    return (
        <div ref={d3Container} />
    );
};

export default NodeDurationGraph;
