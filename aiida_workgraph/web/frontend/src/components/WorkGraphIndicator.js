import React from 'react';
import styled from 'styled-components';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faChevronRight, faCircle } from '@fortawesome/free-solid-svg-icons';

const BreadcrumbContainer = styled.nav`
  background-color: #f7f7f7;
  padding: 10px;
  font-size: 0.9em;
  border-bottom: 1px solid #ddd;
  display: flex;
  align-items: center;
  justify-content: flex-start;
`;

const BreadcrumbLink = styled.a`
  color: #007bff;
  text-decoration: none;
  cursor: pointer;
  margin: 0 5px;
  &:hover {
    text-decoration: underline;
  }
`;

const Separator = styled.span`
  margin: 0 5px;
  color: #ccc;
`;

function Breadcrumbs({ parentWorktrees }) {
  // const separatorIcon = <FontAwesomeIcon icon={faChevronRight} />; // Change this icon as needed
  // const separatorIcon = <FontAwesomeIcon icon={faCircle} />; // Alternative smaller icon
  // const separatorIcon = '>'; // Text-based alternative
  const separatorIcon = '🍞'; // Emoji-based alternative

  return (
    <BreadcrumbContainer aria-label="breadcrumb">
      {parentWorktrees.map(([workgraphName, workgraphId], index) => (
        <React.Fragment key={workgraphId}>
          <BreadcrumbLink href={`/workgraph/${workgraphId}`}>
            {workgraphName}
          </BreadcrumbLink>
          {index < parentWorktrees.length - 1 && <Separator>{separatorIcon}</Separator>}
        </React.Fragment>
      ))}
    </BreadcrumbContainer>
  );
}

export default Breadcrumbs;
