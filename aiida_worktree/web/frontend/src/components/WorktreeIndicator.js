import React from 'react';
import styled from 'styled-components';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faAngleRight } from '@fortawesome/free-solid-svg-icons';

const IndicatorContainer = styled.div`
  background-color: #f7f7f7;
  padding: 10px;
  font-size: 0.9em;
  border-bottom: 1px solid #ddd;
  display: flex;
  align-items: center;
  justify-content: flex-start;
`;

const IndicatorLink = styled.a`
  color: #007bff;
  text-decoration: none;
  cursor: pointer;
  margin-left: 5px;
`;

function WorktreeIndicator({ parentWorktrees }) {
  const hierarchy = parentWorktrees.map(([worktreeName, worktreeId], index) => (
    <span key={worktreeId}>
      {index === 0 && <FontAwesomeIcon size="lg" icon={faAngleRight} />}
      <IndicatorLink href={`/worktree/${worktreeId}`}>{worktreeName}</IndicatorLink>
      {` (${worktreeId})`}
      {index < parentWorktrees.length - 1 && <FontAwesomeIcon size="lg" style={{ marginLeft: '5px' }} icon={faAngleRight} />}
    </span>
  ));

  return <IndicatorContainer>{hierarchy}</IndicatorContainer>;
}

export default WorktreeIndicator;
