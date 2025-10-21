# Project Board Organization

This document describes the recommended project board structure for tracking work in the Bio AI Experiment repository.

## Overview

Project boards help organize and track work across the repository. We recommend using GitHub Projects (new project experience) with the following board structures.

## Recommended Project Boards

### 1. Tool Development Board

**Purpose**: Track the development lifecycle of individual bioinformatics tools

**Columns/Status:**
- 📋 **Backlog** - Tools identified for potential implementation
- 🔍 **Analysis** - Tools being analyzed using TEMPLATE.md
- 📝 **Planning** - Design and planning phase
- 🏗️ **In Progress** - Active development
- 🧪 **Testing** - Testing and validation phase
- 📚 **Documentation** - Documentation being written/updated
- ✅ **Done** - Completed and released
- 🚫 **Won't Do** - Decided not to implement

**Custom Fields:**
- Tool Name (text)
- Category (dropdown): Alignment, QC, Assembly, etc.
- Priority (dropdown): Critical, High, Medium, Low
- Complexity (dropdown): Low, Medium, High
- Status (status)
- Assignee (person)

**Automation:**
- Auto-add issues with "tool-analysis" or "tool-implementation" label
- Move to "In Progress" when issue is assigned
- Move to "Testing" when PR is created
- Move to "Done" when PR is merged

### 2. Feature/Enhancement Board

**Purpose**: Track feature requests and enhancements across all components

**Columns/Status:**
- 📥 **New** - Newly submitted feature requests
- 🔍 **Triage** - Being evaluated for priority and feasibility
- 👍 **Approved** - Accepted for implementation
- 🏗️ **In Progress** - Being implemented
- 👀 **Review** - In code review
- ✅ **Done** - Completed and merged
- ⏸️ **On Hold** - Deferred for later
- 🚫 **Rejected** - Will not implement

**Custom Fields:**
- Component (dropdown): Tool, Library, MCP, Docs, CI
- Priority (dropdown): Critical, High, Medium, Low
- Impact (dropdown): Wide, Moderate, Narrow
- Effort (dropdown): Small, Medium, Large
- Requester (person)
- Status (status)

**Automation:**
- Auto-add issues with "enhancement" label
- Move to "Triage" when created
- Move to "In Progress" when assigned
- Move to "Review" when PR created
- Move to "Done" when PR merged

### 3. Bug Tracking Board

**Purpose**: Track and resolve bugs across the project

**Columns/Status:**
- 🐛 **New** - Newly reported bugs
- 🔍 **Triage** - Evaluating severity and priority
- ✅ **Confirmed** - Bug confirmed and prioritized
- 🏗️ **In Progress** - Being fixed
- 👀 **Review** - Fix in code review
- ✔️ **Resolved** - Fixed and merged
- ❌ **Cannot Reproduce** - Unable to reproduce
- 🔄 **Duplicate** - Duplicate of another issue

**Custom Fields:**
- Severity (dropdown): Critical, High, Medium, Low
- Component (dropdown): Tool name or component
- Affected Version (text)
- Status (status)
- Reporter (person)
- Assignee (person)

**Automation:**
- Auto-add issues with "bug" label
- Move to "Triage" when created
- Move to "In Progress" when assigned
- Move to "Review" when PR created
- Move to "Resolved" when PR merged
- Close issue when moved to "Resolved"

### 4. Research & Analysis Board

**Purpose**: Track analysis work for identifying improvement opportunities

**Columns/Status:**
- 📋 **Proposed** - Tools proposed for analysis
- 🔍 **Analyzing** - Active analysis using TEMPLATE.md
- 📊 **Results** - Analysis completed, decision pending
- ✅ **Approved for Development** - Will be implemented
- 🚫 **Not Recommended** - Will not pursue
- 📚 **Published** - Analysis published and shared

**Custom Fields:**
- Tool Name (text)
- Category (dropdown)
- Priority Score (number)
- Complexity (dropdown)
- Analyst (person)
- Status (status)

**Automation:**
- Auto-add issues with "analysis" or "tool-evaluation" label
- Move to "Analyzing" when assigned
- Track progress through status updates

## Creating Project Boards

### Using GitHub Projects (New Experience)

1. Navigate to the repository
2. Click on "Projects" tab
3. Click "New Project"
4. Select "Board" or "Table" view
5. Name the project (e.g., "Tool Development")
6. Add custom fields as described above
7. Configure automation rules

### Sample Board Views

#### Table View
Provides spreadsheet-like view with all custom fields visible, great for:
- Prioritization and sorting
- Bulk updates
- Filtering by multiple criteria
- Exporting data

#### Board View
Kanban-style board with status columns, great for:
- Visualizing workflow
- Drag-and-drop status updates
- Quick overview of work in progress
- Daily standups

#### Roadmap View
Timeline view showing work over time, great for:
- Release planning
- Dependency tracking
- Long-term planning
- Milestone visualization

## Workflow Integration

### Issue Lifecycle

```
Create Issue → Triage → Planning → Development → Testing → Review → Done
                 ↓          ↓           ↓           ↓        ↓
              Won't Do   On Hold     Blocked      Needs Changes
```

### Automation Rules Examples

1. **Auto-triage**: When issue is created with "bug" label → Move to Bug Tracking Board
2. **Start work**: When issue is assigned → Move to "In Progress"
3. **Review**: When PR is created → Move to "Review"
4. **Complete**: When PR is merged → Move to "Done" and close issue
5. **Blocked**: When "blocked" label added → Move to "Blocked" column

## Project Board Maintenance

### Weekly Tasks
- Review and triage new issues
- Update status of in-progress items
- Clear completed items
- Update priorities based on progress

### Monthly Tasks
- Review "On Hold" items - still relevant?
- Archive old "Done" items
- Update roadmap for next month
- Generate progress reports from board data

## GitHub Project Examples

### Example Queries/Filters

**High Priority Bugs:**
```
is:issue is:open label:bug priority:high,critical
```

**Tools in Development:**
```
is:issue is:open label:tool-implementation status:"In Progress"
```

**This Week's Completions:**
```
is:issue is:closed closed:>=2024-01-15
```

## Integration with Labels

Combine project boards with labels for powerful organization:

- **Type labels**: `bug`, `enhancement`, `documentation`, `analysis`
- **Priority labels**: `priority:critical`, `priority:high`, `priority:medium`, `priority:low`
- **Status labels**: `blocked`, `needs-review`, `good-first-issue`
- **Component labels**: `tool:seqtk`, `tool:prinseq`, `bioformats`, `mcp-server`
- **Area labels**: `performance`, `security`, `testing`, `ci-cd`

## Metrics & Reporting

Track these metrics from project boards:

### Velocity Metrics
- Issues/PRs closed per week
- Average time to close issues
- Average time from "In Progress" to "Done"

### Quality Metrics
- Bug fix rate
- Test coverage increase
- Documentation completeness

### Planning Metrics
- Items in backlog
- Items in progress vs. completed
- Blocked items count

## Team Collaboration

### Daily Standup
Use board views during daily standups:
- What moved forward yesterday?
- What's blocked?
- What's the plan for today?

### Sprint Planning
1. Review completed work from last sprint
2. Evaluate items in backlog
3. Estimate and prioritize next sprint items
4. Assign work to team members
5. Update board with sprint milestone

### Release Planning
1. Filter by milestone
2. Check status of all items
3. Identify blockers
4. Update roadmap view
5. Communicate with stakeholders

## Additional Resources

- [GitHub Projects Documentation](https://docs.github.com/en/issues/planning-and-tracking-with-projects)
- [GitHub Projects Automation](https://docs.github.com/en/issues/planning-and-tracking-with-projects/automating-your-project)
- [Best Practices for Project Management](https://docs.github.com/en/issues/planning-and-tracking-with-projects/learning-about-projects/best-practices-for-projects)

## Getting Started

To set up project boards for this repository:

1. Create the four recommended boards above
2. Configure custom fields for each board
3. Set up automation rules
4. Apply appropriate labels to existing issues
5. Manually add existing issues to relevant boards
6. Train team members on board usage
7. Establish regular board review cadence

---

**Note**: Project boards should be living tools that evolve with the team's needs. Don't hesitate to adjust columns, fields, and workflows as you learn what works best for your team.
