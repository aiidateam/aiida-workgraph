import React from 'react';
import styled from 'styled-components';

const IndicatorContainer = styled.div`
  background-color: #f7f7f7;
  padding: 10px;
  font-size: 0.9em;
  border-bottom: 1px solid #ddd;
  display: flex;
  justify-content: flex-start;
`;

const IndicatorText = styled.span`
  color: #666;
  margin-left: 20px;
`;

const IndicatorLink = styled.a`
  color: #007bff;
  text-decoration: none;
  cursor: pointer;
  margin-left: 20px;
`;

function WorktreeIndicator({ parentWorktrees }) {
  const hierarchy = parentWorktrees.map(([worktreeName, worktreeId], index) => (
    <span key={worktreeId}>
      <IndicatorLink href={`/worktree/${worktreeId}`}>{worktreeName}</IndicatorLink>
      {` (${worktreeId})`}
      {index < parentWorktrees.length - 1 && <IndicatorText> {'->'} </IndicatorText>}
    </span>
  ));

  return <IndicatorContainer>{hierarchy}</IndicatorContainer>;
}

export default WorktreeIndicator;
